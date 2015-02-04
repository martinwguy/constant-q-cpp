/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Constant-Q library
    Copyright (c) 2013-2015 Queen Mary, University of London

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#ifndef CQCHROMAGRAM_H
#define CQCHROMAGRAM_H

#include "CQBase.h"

class CQSpectrogram;

class Chromagram
{
public:
    struct Parameters {
	Parameters(double sr) :
	    sampleRate(sr),
	    lowestOctave(0),
	    octaves(7),
	    bpo(36),
	    tuningFrequency(440.) { }
	double sampleRate;
	int lowestOctave;
	int octaves;
	int bpo;
	double tuningFrequency;
    };

    Chromagram(Parameters params);
    virtual ~Chromagram();

    CQBase::RealBlock process(const CQBase::RealSequence &);
    CQBase::RealBlock getRemainingOutput();

    double getMinFrequency() const { return m_minFrequency; }
    double getMaxFrequency() const { return m_maxFrequency; }

    std::string getBinName(int bin) const;
    
    bool isValid() const;
    int getColumnHop() const;
    int getLatency() const;
    
private:
    Parameters m_params;
    CQSpectrogram *m_cq;
    double m_minFrequency;
    double m_maxFrequency;
    CQBase::RealBlock convert(const CQBase::RealBlock &);
};

#endif


    
