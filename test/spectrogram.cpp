
// #include "ConstantQ.h"
#include "CQSpectrogram.h"

#include <sndfile.h>
#include <png.h>

#include <iostream>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

#include <cstring>

#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <cstdlib>

static void output_frames(CQSpectrogram::RealBlock rb);
static void write_png(char *fileNameOut);

// Number of data items in each analysis frame
static unsigned int nbuckets;
// Total number of frames to write to output file
static unsigned int nframes;

// How many of the initial output frames do we still need to ignore
// to compensate for the latency?
static png_uint_32 framesToSkip;
// How many more frames do we need to output into the PNG file?
static png_uint_32 framesToOutput;

/* Data is generated column-by-column but PNGs need writing
 * row-by-row so store output columns as they are generated
 * then output them the right way up at the end.
 * This also allows us to find the maximum amplitude and
 * normalize the output so that the loudest pixel is white..
 */
static float **frames;	// All the output frames
static float maxAmp = -1.0/0.0;	/* -infinity */
static float dynRange = 100;	/* dynamic range of output in decibels */

int main(int argc, char **argv)
{
    double maxFreq = 0;
    double minFreq = 0;
    double atomHopFactor = 0.0; /* If 0, not given as a parameter */
    int bpo = 0;
    bool help = false;
    CQSpectrogram::Interpolation interpolation = CQSpectrogram::InterpolateLinear;
    CQParameters::WindowType window = CQParameters::Hann;

    // How many output frames have process() and getRemaining() returned to us?
    int framesReceived = 0;

    int c;

    while (1) {
	int optionIndex = 0;
	
	static struct option longOpts[] = {
	    { "help", 0, 0, 'h', },
	    { "maxfreq", 1, 0, 'x', },
	    { "minfreq", 1, 0, 'n', },
	    { "dynrange", 1, 0, 'd', },
	    { "atomhopfactor", 1, 0, 'a', },
	    { "bpo", 1, 0, 'b' },
	    { "interpolation", 1, 0, 'i' },
	    { "window", 1, 0, 'w' },
	    { 0, 0, 0, 0 },
	};

	c = getopt_long(argc, argv,
			"hx:n:d:a:b:i:w:",
			longOpts, &optionIndex);
	if (c == -1) break;

	switch (c) {
	case 'h': help = true; break;
	case 'x': maxFreq = atof(optarg); break;
	case 'n': minFreq = atof(optarg); break;
	case 'd': dynRange = atof(optarg); break;
	case 'a': atomHopFactor = atof(optarg); break;
	case 'b': bpo = atoi(optarg); break;
	case 'i': switch (optarg[0]) {
		  case 'h': interpolation = CQSpectrogram::InterpolateHold; break;
		  case 'l': interpolation = CQSpectrogram::InterpolateLinear; break;
		  case 'z': interpolation = CQSpectrogram::InterpolateZeros; break;
		  default: help = true; break;
		  } break;
	case 'w': switch (optarg[0]) {
		  case 'n': window = CQParameters::Hann; break;
		  case 'k': window = CQParameters::Blackman; break;
		  case 'h': window = CQParameters::BlackmanHarris; break;
		  case 'N': window = CQParameters::SqrtHann; break;
		  case 'K': window = CQParameters::SqrtBlackman; break;
		  case 'H': window = CQParameters::SqrtBlackmanHarris; break;
		  default: help = true; break;
		  } break;
	default: help = true; break;
	}
    }

    if (help || (optind + 2 != argc && optind + 3 != argc)) {
	cerr << endl;
	cerr << "Usage: " << argv[0] << " [options] infile.wav outfile.png" << endl;
	cerr << endl;
	cerr << "Options:" << endl;
	cerr << "  -x<X>, --maxfreq <X>  Maximum frequency (default = sample rate / 3)" << endl;
	cerr << "  -n<X>, --minfreq <X>  Minimum frequency (default = 100, actual min may vary)" << endl;
	cerr << "  -d<X>, --dynrange <X> Dynamic range of output in dB (default: 100)" << endl;
	cerr << "  -a<X>, --atomhopfactor <X> Atom hop factor (default = 0.25)" << endl;
	cerr << "  -b<X>, --bpo <X>      Bins per octave   (default = 60)" << endl;
	cerr << "  -i<X>, --interpolation <X>" << endl;
	cerr << "                        Interpolation type: h = Hold, l = Linear, Z = Zeros" << endl;
	cerr << "  -w<X>, --window <X>   Window type n = Hann, k = Blackman, h = BlackmanHarris" << endl;
	cerr << "                        Use N, K and H for the square root of the above." << endl;
	cerr << "  -h, --help            Print this help" << endl;
	cerr << endl;
	cerr << "This program performs a Constant-Q spectrogram with the" << endl;
	cerr << "requested parameters and writes the results to a PNG file." << endl;
	cerr << endl;
	return 2;
    }

    char *fileName = strdup(argv[optind++]);
    char *fileNameOut = strdup(argv[optind++]);

    SNDFILE *sndfile;
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    /* Open sound file for reading */
    sndfile = sf_open(fileName, SFM_READ, &sfinfo);
    if (!sndfile) {
	cerr << "ERROR: Failed to open input file \"" << fileName << "\": "
	     << sf_strerror(sndfile) << endl;
	return 1;
    }

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (minFreq == 0.0) minFreq = 100;
    if (bpo == 0) bpo = 60;

    CQParameters params(sfinfo.samplerate, minFreq, maxFreq, bpo);
    params.window = window;
    if (atomHopFactor != 0.0) params.atomHopFactor = atomHopFactor;
    CQSpectrogram cq(params, interpolation);

    nbuckets = (png_uint_32) (log2 (maxFreq / minFreq) * bpo + 0.5);
    nframes = sfinfo.frames / cq.getColumnHop();

    frames = new float * [nframes];
    if (!frames) {
	cerr << "Out of memory." << endl;
	exit(1);
    }

    int inframe = 0;

    vector<double> buffer;

    framesToSkip = (png_uint_32)((double)cq.getLatency() / cq.getColumnHop() + 0.5);
    framesToOutput = nframes;

    while (inframe < sfinfo.frames) {
        int count = -1;
        CQSpectrogram::RealBlock rb;
	
	if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0) {
	    break;
	}

	// Convert multi-channel sound files to mono
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

	rb = cq.process(cqin);
	// RealBlock is a vector of RealColumn's and
	// RealColumn is a vector of double.
        if (rb.size() > 0) {
	    framesReceived += rb.size();
	    output_frames(rb);
	}

	inframe += count;
    }

    sf_close(sndfile);

    CQSpectrogram::RealBlock rb = cq.getRemainingOutput();
    if (rb.size() > 0) {
	framesReceived += rb.size();
	output_frames(rb);
    }

    write_png(fileNameOut);

    return 0;
}

static void
output_frames(CQSpectrogram::RealBlock rb)
{
    static float **framesp = NULL;
    if (framesp == NULL) framesp = frames;

    for (unsigned i = 0; i < rb.size(); i++) {
	if (framesToSkip > 0) {
	    framesToSkip--;
	    continue;
        }
        if (framesToOutput <= 0)
	    return;

        float *frame = new float[nbuckets];
        if (!frame) {
	    cerr << "Out of memory." << endl;
	    exit(1);
        }
	for (unsigned j = 0; j < nbuckets; j++) {
	    float value = 20.0 * log10(rb[i][j]);	// Convert to dB
	    frame[j] = value;
	    if (value > maxAmp) maxAmp=value;
	}
	*framesp++ = frame;

        framesToOutput--;
    }
}

static void
write_png(char *fileNameOut)
{
    //
    // Write PNG header
    //
    png_structp png_ptr;
    png_infop png_info;
    png_uint_32 png_width = nframes;
    png_uint_32 png_height = nbuckets;

    /* Create PNG file for writing */
    FILE *png_fp = fopen(fileNameOut, "wb");
    if (!png_fp) {
	cerr << "Cannot create " << fileNameOut << endl;
	exit(1);
    }
    png_ptr = png_create_write_struct(
	PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
	cerr << "Cannot create PNG writer" << endl;
	exit(1);
    }
    png_info = png_create_info_struct(png_ptr);
    if (!png_info) {
	cerr << "Cannot create PNG info" << endl;
	exit(1);
    }
    if (setjmp(png_jmpbuf(png_ptr)))
    {
	cerr << "Something went wrong in the PNG-writing subsystem. Quitting..." << endl;
	png_destroy_write_struct(&png_ptr, &png_info);
        fclose(png_fp);
	exit(1);
    }
    png_init_io(png_ptr, png_fp);
    png_set_IHDR(png_ptr, png_info,
	png_width, png_height,
	8, // bit depth
	PNG_COLOR_TYPE_GRAY,
	PNG_INTERLACE_NONE,
	PNG_COMPRESSION_TYPE_DEFAULT,
	PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, png_info);

    //
    // Write image
    //
    png_bytep png_row = new png_byte[png_width];
    if (!png_row) {
	cerr << "Out of memory." << endl;
	exit(1);
    }

    for (png_uint_32 y=0; y < png_height; y++) {
        for (png_uint_32 x=0; x < png_width; x++) {
	    float value = (maxAmp - frames[x][y]) / dynRange;
	    if (value < 0.0) { value = 0.0; }
	    if (value > 1.0) { value = 1.0; }
	    png_row[x] = ((1.0-value) * 255) + 0.5;
	}
	png_write_row(png_ptr, png_row);
    }

    png_write_end(png_ptr, NULL);
}
