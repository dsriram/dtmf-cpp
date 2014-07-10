//
// Utilize the DtmfDetector to detect tones in an AU file.
// The file must be 8KHz, PCM encoded, mono.
//

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdint.h>

#include "DtmfDetector.hpp"

//
// The string ".snd" in big-endian byte ordering.  This identifies the file as
// an AU sound file.
//
#define AU_MAGIC 0x2e736e64

//
// The size of the buffer we use for reading & processing the audio samples.
//
#define BUFLEN 512

using namespace std;

//
// Swap the endianness of an 4-byte integer.
//

//////////////////////////////////////////////////////////////
//  Filter Code Definitions
//////////////////////////////////////////////////////////////
 
// maximum number of inputs that can be handled
// in one function call
#define MAX_INPUT_LEN   512
// maximum length of filter than can be handled
#define MAX_FLT_LEN     63
// buffer to hold all of the input samples
#define BUFFER_LEN      (MAX_FLT_LEN - 1 + MAX_INPUT_LEN)
 
// array to hold input samples
double insamp[ BUFFER_LEN ];
 
// FIR init
void firFloatInit( void )
{
    memset( insamp, 0, sizeof( insamp ) );
}
 
// the FIR filter function
void firFloat( double *coeffs, double *input, double *output,
       int length, int filterLength )
{
    double acc;     // accumulator for MACs
    double *coeffp; // pointer to coefficients
    double *inputp; // pointer to input samples
    int n;
    int k;
 
    // put the new samples at the high end of the buffer
    memcpy( &insamp[filterLength - 1], input,
            length * sizeof(double) );
    // apply the filter to each input sample
    for ( n = 0; n < length; n++ ) {
        // calculate output n
        coeffp = coeffs;
        inputp = &insamp[filterLength - 1 + n];
        acc = 0;
        for ( k = 0; k < filterLength; k++ ) {
            acc += (*coeffp++) * (*inputp--);
        }
        output[n] = acc;
    }
    // shift input samples back in time for next time
    memmove( &insamp[0], &insamp[length],
            (filterLength - 1) * sizeof(double) );
 
}
 
//////////////////////////////////////////////////////////////
//  Test program
//////////////////////////////////////////////////////////////
 
// bandpass filter centred around 1000 Hz
// sampling rate = 8000 Hz
 
#define FILTER_LEN  63
double coeffs[ FILTER_LEN ] =
{
  -0.000024, 
-0.000002, 
-0.000034, 
-0.000279, 
-0.000476, 
-0.000022, 
0.000923, 
0.001128, 
0.000261, 
0.000274, 
0.002470, 
0.003720, 
-0.000312, 
-0.007125, 
-0.008251, 
-0.002189, 
0.000169, 
-0.007995, 
-0.013064, 
0.003114, 
0.030448, 
0.035229, 
0.009407, 
-0.007612, 
0.014938, 
0.035829, 
-0.017542, 
-0.131490, 
-0.178252, 
-0.056293, 
0.159365, 
0.267500, 
0.159365, 
-0.056293, 
-0.178252, 
-0.131490, 
-0.017542, 
0.035829, 
0.014938, 
-0.007612, 
0.009407, 
0.035229, 
0.030448, 
0.003114, 
-0.013064, 
-0.007995, 
0.000169, 
-0.002189, 
-0.008251, 
-0.007125, 
-0.000312, 
0.003720, 
0.002470, 
0.000274, 
0.000261, 
0.001128, 
0.000923, 
-0.000022, 
-0.000476, 
-0.000279, 
-0.000034, 
-0.000002, 
-0.000024
};
 
void intToFloat( int16_t *input, double *output, int length )
{
    int i;
 
    for ( i = 0; i < length; i++ ) {
        output[i] = (double)input[i];
    }
}
 
void floatToInt( double *input, int16_t *output, int length )
{
    int i;
 
    for ( i = 0; i < length; i++ ) {
        if ( input[i] > 32767.0 ) {
            input[i] = 32767.0;
        } else if ( input[i] < -32768.0 ) {
            input[i] = -32768.0;
        }
        // convert
        output[i] = (int16_t)input[i];
    }
}



uint32_t
swap32(uint32_t a)
{
    char b1, b2, b3, b4 = 0;
    b1 = a;
    b2 = a >> 8;
    b3 = a >> 16;
    b4 = a >> 24;
    return (b1 << 24) + (b2 << 16) + (b3 << 8) + b4;
}

struct au_header
{
    uint32_t magic;
    uint32_t header_size;
    uint32_t nsamples;
    uint32_t encoding;
    uint32_t sample_rate;
    uint32_t nchannels;
};

string
au_header_tostr(au_header &h)
{
    stringstream ss;
    ss << h.header_size << " header bytes, " << 
        h.nsamples << " samples, encoding type: " << h.encoding << ", " 
        << h.sample_rate << "Hz, " << h.nchannels << " channels";
    return ss.str();
}

int
main(int argc, char **argv)
{
    if (argc != 2)
    {
        cerr << "usage: " << argv[0] << " filename.au" << endl;
        return 1;
    }

    ifstream fin(argv[1], ios::binary);
    if (!fin.good())
    {
        cerr << argv[1] << ": unable to open file" << endl;
        return 1;
    }
    au_header header;
    fin.read((char *)&header, sizeof(au_header));

    //
    // The data in the AU file header is stored in big-endian byte ordering.
    // If we're on a little-endian machine, we need to reorder the bytes.
    // While other arrangements (e.g. middle-endian) are also technically
    // possible, this example does not support them, since they are relatively
    // rare.
    //
    if (header.magic == AU_MAGIC)
    {
    }
    else if (header.magic == swap32(AU_MAGIC))
    {
        header.header_size = swap32(header.header_size);
        header.nsamples = swap32(header.nsamples);
        header.encoding = swap32(header.encoding);
        header.sample_rate = swap32(header.sample_rate);
        header.nchannels = swap32(header.nchannels);
    }
    else
    {
        cerr << "bad magic number: " << hex << header.magic << endl;
        return 1;
    }

    cout << argv[1] << ": " << au_header_tostr(header) << endl;
    //
    // This example only supports a specific type of AU format:
    //
    // - no additional data in the header
    // - linear PCM encoding
    // - 8KHz sample rate
    // - mono
    //
    if 
    (
        header.header_size != 24
        ||
        header.encoding != 2
        ||
        header.sample_rate != 8000
        ||
        header.nchannels != 1
    )
    {
        cerr << argv[1] << ": unsupported AU format" << endl;
        return 1;
    }

    short input[BUFLEN];
    short output[BUFLEN];
    double floatInput[BUFLEN];
    double floatOutput[BUFLEN];
    char cbuf[BUFLEN];
    short sbuf[BUFLEN];
    DtmfDetector detector(BUFLEN);
    firFloatInit();
    for (uint32_t i = 0; i < header.nsamples; i += BUFLEN)
    {
        fin.read(cbuf, BUFLEN);
        //
        // Promote our 8-bit samples to 16 bits, since that's what the detector
        // expects.  Shift them left during promotion, since the decoder won't
        // pick them up otherwise (volume too low).
        //
        for (int j = 0; j < BUFLEN; ++j)
            sbuf[j] = cbuf[j] << 8;
        detector.zerosIndexDialButton();

        intToFloat( sbuf, floatInput, BUFLEN );
        // perform the filtering
        firFloat( coeffs, floatInput, floatOutput, BUFLEN,
               FILTER_LEN );
        // convert to ints
        floatToInt( floatOutput, sbuf, BUFLEN );

         //detector.dtmfDetecting(sbuf);
         cout << i << ": `" << detector.DTMF_detection(sbuf) << "'" << endl;
         //cout << i << ": `" << detector.getDialButtonsArray() << "'" << endl;
    }
    cout << endl;
    fin.close();

    return 0;
}
