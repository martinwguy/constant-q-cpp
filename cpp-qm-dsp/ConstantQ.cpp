
#include "ConstantQ.h"

#include "CQKernel.h"

#include "dsp/rateconversion/Resampler.h"
#include "maths/MathUtilities.h"
#include "dsp/transforms/FFT.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::complex;
using std::cerr;
using std::endl;

typedef std::complex<double> C;

ConstantQ::ConstantQ(double sampleRate,
                     double minFreq,
                     double maxFreq,
                     int binsPerOctave) :
    m_sampleRate(sampleRate),
    m_maxFrequency(maxFreq),
    m_minFrequency(minFreq),
    m_binsPerOctave(binsPerOctave),
    m_fft(0)
{
    if (minFreq <= 0.0 || maxFreq <= 0.0) {
        throw std::invalid_argument("Frequency extents must be positive");
    }

    initialise();
}

ConstantQ::~ConstantQ()
{
    delete m_fft;
    for (int i = 0; i < (int)m_decimators.size(); ++i) {
        delete m_decimators[i];
    }
    delete m_kernel;
}

double
ConstantQ::getMinFrequency() const
{
    return m_p.minFrequency / pow(2.0, m_octaves - 1);
}

double
ConstantQ::getBinFrequency(int bin) const
{
    return getMinFrequency() * pow(2, (double(bin) / getBinsPerOctave()));
}

void
ConstantQ::initialise()
{
    m_octaves = int(ceil(log2(m_maxFrequency / m_minFrequency)));
    m_kernel = new CQKernel(m_sampleRate, m_maxFrequency, m_binsPerOctave);
    m_p = m_kernel->getProperties();
    
    // Use exact powers of two for resampling rates. They don't have
    // to be related to our actual samplerate: the resampler only
    // cares about the ratio, but it only accepts integer source and
    // target rates, and if we start from the actual samplerate we
    // risk getting non-integer rates for lower octaves

    int sourceRate = pow(2, m_octaves);
    vector<int> latencies;

    // top octave, no resampling
    latencies.push_back(0);
    m_decimators.push_back(0);

    for (int i = 1; i < m_octaves; ++i) {

        int factor = pow(2, i);

        Resampler *r = new Resampler
            (sourceRate, sourceRate / factor, 60, 0.02);

        // We need to adapt the latencies so as to get the first input
        // sample to be aligned, in time, at the decimator output
        // across all octaves.
        // 
        // Our decimator uses a linear phase filter, but being causal
        // it is not zero phase: it has a latency that depends on the
        // decimation factor. Those latencies have been calculated
        // per-octave and are available to us in the latencies
        // array. Left to its own devices, the first input sample will
        // appear at output sample 0 in the highest octave (where no
        // decimation is needed), sample number latencies[1] in the
        // next octave down, latencies[2] in the next one, etc. We get
        // to apply some artificial per-octave latency after the
        // decimator in the processing chain, in order to compensate
        // for the differing latencies associated with different
        // decimation factors. How much should we insert?
        //
        // The outputs of the decimators are at different rates (in
        // terms of the relation between clock time and samples) and
        // we want them aligned in terms of time. So, for example, a
        // latency of 10 samples with a decimation factor of 2 is
        // equivalent to a latency of 20 with no decimation -- they
        // both result in the first output sample happening at the
        // same equivalent time in milliseconds.
	// 
	// So here we record the latency added by the decimator, in
	// terms of the sample rate of the undecimated signal. Then we
	// use that to compensate in a moment, when we've discovered
	// what the longest latency across all octaves is.

        latencies.push_back(r->getLatency() * factor);
        m_decimators.push_back(r);
    }

    m_bigBlockSize = m_p.fftSize * pow(2, m_octaves - 1);

    // Now add in the extra padding and compensate for hops that must
    // be dropped in order to align the atom centres across
    // octaves. Again this is a bit trickier because we are doing it
    // at input rather than output and so must work in per-octave
    // sample rates rather than output blocks

    int emptyHops = m_p.firstCentre / m_p.atomSpacing;

    vector<int> drops;
    for (int i = 0; i < m_octaves; ++i) {
	int factor = pow(2, i);
	int dropHops = emptyHops * pow(2, m_octaves - i - 1) - emptyHops;
	int drop = ((dropHops * m_p.fftHop) * factor) / m_p.atomsPerFrame;
	drops.push_back(drop);
    }

    int maxLatPlusDrop = 0;
    for (int i = 0; i < m_octaves; ++i) {
	int latPlusDrop = latencies[i] + drops[i];
	if (latPlusDrop > maxLatPlusDrop) maxLatPlusDrop = latPlusDrop;
    }

    // we want to design totalLatency such that totalLatency -
    // latencies[0] - drops[0] is a multiple of m_p.fftHop, so that we
    // can get identical results in octave 0 to our reference
    // implementation, making for easier testing (though other octaves
    // will differ because of different resampler implementations)

    int totalLatency = maxLatPlusDrop;
    int lat0 = totalLatency - latencies[0] - drops[0];
    totalLatency = ceil(double(lat0 / m_p.fftHop) * m_p.fftHop)
	+ latencies[0] + drops[0];

//    cerr << "total latency = " << totalLatency << endl;

    // Padding as in the reference (will be introduced with the
    // latency compensation in the loop below)
    m_outputLatency = totalLatency + m_bigBlockSize
	- m_p.firstCentre * pow(2, m_octaves-1);

//    cerr << "m_bigBlockSize = " << m_bigBlockSize << ", firstCentre = "
//	 << m_p.firstCentre << ", m_octaves = " << m_octaves << ", so m_outputLatency = " << m_outputLatency << endl;

    for (int i = 0; i < m_octaves; ++i) {

	double factor = pow(2, i);

	// Calculate the difference between the total latency applied
	// across all octaves, and the existing latency due to the
	// decimator for this octave, and then convert it back into
	// the sample rate appropriate for the output latency of this
	// decimator -- including one additional big block of padding
	// (as in the reference).

	double octaveLatency =
	    double(totalLatency - latencies[i] - drops[i]
		   + m_bigBlockSize) / factor;

        m_buffers.push_back
            (vector<double>(int(round(octaveLatency)), 0.0));
    }

    m_fft = new FFTReal(m_p.fftSize);
}

vector<vector<double> > 
ConstantQ::process(const vector<double> &td)
{
    m_buffers[0].insert(m_buffers[0].end(), td.begin(), td.end());

    for (int i = 1; i < m_octaves; ++i) {
        vector<double> dec = m_decimators[i]->process(td.data(), td.size());
        m_buffers[i].insert(m_buffers[i].end(), dec.begin(), dec.end());
    }

    vector<vector<double> > out;

    while (true) {

	// We could have quite different remaining sample counts in
	// different octaves, because (apart from the predictable
	// added counts for decimator output on each block) we also
	// have variable additional latency per octave
	bool enough = true;
	for (int i = 0; i < m_octaves; ++i) {
	    int required = m_p.fftSize * pow(2, m_octaves - i - 1);
//	    cerr << "for octave " << i << ", buf len =  "<< m_buffers[i].size() << " (need " << required << ")" << endl;

	    if ((int)m_buffers[i].size() < required) {
		enough = false;
	    }
	}
	if (!enough) break;

        int base = out.size();
        int totalColumns = pow(2, m_octaves - 1) * m_p.atomsPerFrame;
        for (int i = 0; i < totalColumns; ++i) {
            out.push_back(vector<double>(m_p.binsPerOctave * m_octaves, 0.0));
        }

        for (int octave = 0; octave < m_octaves; ++octave) {

            int blocksThisOctave = pow(2, (m_octaves - octave - 1));

            for (int b = 0; b < blocksThisOctave; ++b) {
                vector<vector<double> > block = processOctaveBlock(octave);
                
                for (int j = 0; j < m_p.atomsPerFrame; ++j) {

		    for (int k = 0; k < pow(2, octave); ++k) {

			int target = base + k +
			    (b * (totalColumns / blocksThisOctave) + 
			     (j * ((totalColumns / blocksThisOctave) /
				   m_p.atomsPerFrame)));

			for (int i = 0; i < m_p.binsPerOctave; ++i) {
			    out[target][m_p.binsPerOctave * octave + i] = 
				block[j][m_p.binsPerOctave - i - 1];
			}
		    }
                }
            }
        }
    }

    return out;
}

vector<vector<double> >
ConstantQ::getRemainingBlocks()
{
    // Same as padding added at start, though rounded up
    int pad = ceil(double(m_outputLatency) / m_bigBlockSize) * m_bigBlockSize;
    vector<double> zeros(pad, 0.0);
    return process(zeros);
}

vector<vector<double> >
ConstantQ::processOctaveBlock(int octave)
{
    vector<double> ro(m_p.fftSize, 0.0);
    vector<double> io(m_p.fftSize, 0.0);

    m_fft->forward(m_buffers[octave].data(), ro.data(), io.data());

    vector<double> shifted;
    shifted.insert(shifted.end(), 
                   m_buffers[octave].begin() + m_p.fftHop,
                   m_buffers[octave].end());
    m_buffers[octave] = shifted;

    vector<C> cv;
    for (int i = 0; i < m_p.fftSize; ++i) {
        cv.push_back(C(ro[i], io[i]));
    }

    vector<C> cqrowvec = m_kernel->process(cv);

    // Reform into a column matrix
    vector<vector<double> > cqblock;
    for (int j = 0; j < m_p.atomsPerFrame; ++j) {
        cqblock.push_back(vector<double>());
        for (int i = 0; i < m_p.binsPerOctave; ++i) {
            cqblock[j].push_back(abs(cqrowvec[i * m_p.atomsPerFrame + j]));
        }
    }

    return cqblock;
}


