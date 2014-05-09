
#include "ConstantQ.h"
#include "CQInverse.h"

#include <sndfile.h>

#include <iostream>

using std::vector;
using std::cerr;
using std::endl;

#include <cstring>

#include <unistd.h>
#include <sys/time.h>

int main(int argc, char **argv)
{
    if (argc != 3) {
	cerr << "Usage: " << argv[0] << " infile.wav outfile.wav" << endl;
	return 2;
    }

    char *fileName = strdup(argv[1]);
    char *fileNameOut = strdup(argv[2]);

    SNDFILE *sndfile;
    SNDFILE *sndfileOut;
    SF_INFO sfinfo;
    SF_INFO sfinfoOut;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(fileName, SFM_READ, &sfinfo);
    if (!sndfile) {
	cerr << "ERROR: Failed to open input file \"" << fileName << "\": "
	     << sf_strerror(sndfile) << endl;
	return 1;
    }

    sfinfoOut.channels = 1;
    sfinfoOut.format = sfinfo.format;
    sfinfoOut.frames = sfinfo.frames;
    sfinfoOut.samplerate = sfinfo.samplerate;
    sfinfoOut.sections = sfinfo.sections;
    sfinfoOut.seekable = sfinfo.seekable;

    sndfileOut = sf_open(fileNameOut, SFM_WRITE, &sfinfoOut) ;
    if (!sndfileOut) {
	cerr << "ERROR: Failed to open output file \"" << fileNameOut << "\" for writing: "
	     << sf_strerror(sndfileOut) << endl;
	return 1;
    }

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    ConstantQ cq(sfinfo.samplerate,
		 100, sfinfo.samplerate / 3,
		 60);

    CQInverse cqi(sfinfo.samplerate,
		 100, sfinfo.samplerate / 3,
		 60);

    int inframe = 0;
    int outframe = 0;
    int latency = cq.getLatency() + cqi.getLatency();

    cerr << "forward latency = " << cq.getLatency() << ", inverse latency = " 
	 << cqi.getLatency() << ", total = " << latency << endl;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    while (inframe < sfinfo.frames) {

        int count = -1;
	
	if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0) {
	    break;
	}

	vector<double> cqin;
	for (int i = 0; i < count; ++i) {
	    double v = fbuf[i * channels];
	    if (channels > 1) {
		for (int c = 1; c < channels; ++c) {
		    v += fbuf[i * channels + c];
		}
		v /= channels;
	    }
	    cqin.push_back(v);
	}
	
	vector<double> cqout = cqi.process(cq.process(cqin));

	if (outframe >= latency) {

	    sf_writef_double(sndfileOut, 
			     cqout.data(), 
			     cqout.size());

	} else if (outframe + (int)cqout.size() >= latency) {

	    int offset = latency - outframe;
	    sf_writef_double(sndfileOut, 
			     cqout.data() + offset,
			     cqout.size() - offset);
	}

	inframe += count;
	outframe += cqout.size();
    }

    vector<double> r1 = cqi.process(cq.getRemainingOutput());
    vector<double> r2 = cqi.getRemainingOutput();

    sf_writef_double(sndfileOut, r1.data(), r1.size());
    sf_writef_double(sndfileOut, r2.data(), r2.size());

    outframe += r1.size();
    outframe += r2.size();

    sf_close(sndfile);
    sf_close(sndfileOut);

    cerr << "in: " << inframe << ", out: " << outframe - latency << endl;

    timeval etv;
    (void)gettimeofday(&etv, 0);
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
	etv.tv_usec += 1000000;
	etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << endl;

    return 0;
}

