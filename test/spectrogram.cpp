
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

static void output_rows(png_structp png_ptr,
			png_uint_32 png_width,
			CQSpectrogram::RealBlock rb,
			double maxAmp);

static void show_min_max(void);

// How many of the initial output columns do we still need to ignore
// to compensate for the latency?
static png_uint_32 columnsToSkip;
// How many more columns do we need to output into the PNG file?
static png_uint_32 columnsToOutput;

int main(int argc, char **argv)
{
    double maxFreq = 0;
    double minFreq = 0;
    double maxAmp = 1.0;
    double atomHopFactor = 0.0; /* If 0, not given as a parameter */
    int bpo = 0;
    bool help = false;
    CQSpectrogram::Interpolation interpolation = CQSpectrogram::InterpolateLinear;
    CQParameters::WindowType window = CQParameters::Hann;

    // How many output columns have process() and getRemaining() returned to us?
    int columnsReceived = 0;

    int c;

    while (1) {
	int optionIndex = 0;
	
	static struct option longOpts[] = {
	    { "help", 0, 0, 'h', },
	    { "maxfreq", 1, 0, 'x', },
	    { "minfreq", 1, 0, 'n', },
	    { "maxamp", 1, 0, 'm', },
	    { "atomhopfactor", 1, 0, 'a', },
	    { "bpo", 1, 0, 'b' },
	    { "interpolation", 1, 0, 'i' },
	    { "window", 1, 0, 'w' },
	    { 0, 0, 0, 0 },
	};

	c = getopt_long(argc, argv,
			"hx:n:m:a:b:i:w:",
			longOpts, &optionIndex);
	if (c == -1) break;

	switch (c) {
	case 'h': help = true; break;
	case 'x': maxFreq = atof(optarg); break;
	case 'n': minFreq = atof(optarg); break;
	case 'm': maxAmp = atof(optarg); break;
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
	cerr << "  -m<X>, --maxamp <X>   Predited maximum output value (default = 1.0)" << endl;
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

#if 0
    cerr << "octave boundaries: ";
    for (int i = 0; i < cq.getOctaves(); ++i) {
        cerr << cq.getMaxFrequency() / pow(2, i) << " ";
    }
    cerr << endl;
#endif

    png_structp png_ptr;
    png_infop png_info;
    png_uint_32 png_width;

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
    png_width = (png_uint_32) (log2 (maxFreq / minFreq) * bpo + 0.5);
    png_uint_32 png_height = sfinfo.frames / cq.getColumnHop();
    png_init_io(png_ptr, png_fp);
    png_set_IHDR(png_ptr, png_info,
	png_width, png_height,
	8, // bit depth
	PNG_COLOR_TYPE_GRAY,
	PNG_INTERLACE_NONE,
	PNG_COMPRESSION_TYPE_DEFAULT,
	PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, png_info);

    int inframe = 0;

    vector<double> buffer;

    columnsToSkip = (png_uint_32)((double)cq.getLatency() / cq.getColumnHop() + 0.5);
    columnsToOutput = png_height;

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
	    columnsReceived += rb.size();
	    output_rows(png_ptr, png_width, rb, maxAmp);
	}

	inframe += count;
    }

    sf_close(sndfile);

    CQSpectrogram::RealBlock rb = cq.getRemainingOutput();
    if (rb.size() > 0) {
	columnsReceived += rb.size();
	output_rows(png_ptr, png_width, rb, maxAmp);
    }

    png_write_end(png_ptr, NULL);

    show_min_max();

    return 0;
}

static double minval, maxval;
static bool initval=0;

static void
show_min_max()
{
    cerr << "Maximum output value: " << maxval << endl;
}

static void
output_rows(png_structp png_ptr, png_uint_32 png_width, CQSpectrogram::RealBlock rb, double maxAmp)
{
    static png_bytep png_row = NULL;
    if (png_row == NULL) png_row = (png_bytep) malloc(png_width);
    if (!png_row) {
	cerr << "Out of memory." << endl;
	exit(1);
    }

    if (!initval) minval = maxval = rb[0][0];

    /* Convert doubles to 0-255 */
    for (unsigned i = 0; i < rb.size(); i++) {
	if (columnsToSkip > 0) {
	    columnsToSkip--;
	    continue;
        }
        if (columnsToOutput <= 0)
	    return;
	for (unsigned from = 0, to = png_width; to > 0; ) {
	    double value = rb[i][from++];
	    if (value < minval) minval=value;
	    if (value > maxval) maxval=value;
	    if (value < 0.0) { value = 0.0; }
	    if (value > maxAmp ) { value = maxAmp; }
	    png_row[--to] = (value / maxAmp * 255) + 0.5;
	}
	png_write_row(png_ptr, png_row);
	columnsToOutput--;
    }
}
