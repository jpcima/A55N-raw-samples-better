#!/usr/bin/env dub
/+ dub.sdl: name "hello" ; libs "sndfile" +/

import sndfile;
import std.stdio;
import std.format;
import std.string;
import std.regex;
import std.conv;
import std.path;
import std.math;
import std.getopt;
import std.algorithm;
import std.exception;

struct Sound
{
    ulong frames;
    float samplerate;
    float[][] channels;
};

struct MidiRange
{
    uint basekey;
    uint numkeys;
};

void main(string[] args)
{
    string inputPath;
    string outputPath;

    auto helpInfo = getopt(args,
        "i|input", &inputPath,
        "o|output", &outputPath);

    if (helpInfo.helpWanted)
    {
        defaultGetoptPrinter("Some information about the program.", helpInfo.options);
        return;
    }

    if (args.length != 1)
    {
        throw new Exception("Invalid arguments");
    }
    if (inputPath.empty)
    {
        throw new Exception("No input file given");
    }
    if (outputPath.empty)
    {
        throw new Exception("No output directory given");
    }

    MidiRange midiRange = getMidiKeysOfFile(inputPath);

    Sound *snd = readSoundFile(inputPath);
    float samplerate = snd.samplerate;

    float[][] seqs = extractNonZeroSequences(snd.channels[0], samplerate, 0.1);

    if (seqs.length != midiRange.numkeys)
        throw new Exception("Bad number of segments");

    seqs = reduceToCenterSamples(seqs, samplerate, 2.0);

    ///
    if (false)
    {
        ulong j = 0;
        for (ulong c = 0; c < seqs.length; ++c)
        {
            float[] seq = seqs[c];
            bool[] peaks = new bool[seq.length];
            foreach (ulong peak; findPeaks(seq))
            {
                peaks[peak] = true;
            }
            for (ulong i = 0; i < seq.length; ++i)
            {
                writefln!"%e %e %e"(j / samplerate, seq[i], peaks[i] ? seq[i] : 0.0);
                ++j;
            }
        }
    }

    /// extract the periodic part of the signal
    for (ulong iseq = 0; iseq < seqs.length; ++iseq)
    {
        uint currentNote = midiRange.basekey + cast(uint)iseq;
        float currentFreq = 440.0f * pow(2.0f, (currentNote - 69.0f) / 12.0f);

        stderr.writefln!"";
        stderr.writefln!"    Wave %d key=%d freq=%f"(iseq, currentNote, currentFreq);

        //------------------------------

        float[] nonPeriodic = seqs[iseq];

        ulong[] pks = findPeaks(nonPeriodic);
        if (pks.empty)
            throw new Exception("Wave does not have any peaks");

        // find highest peak
        ulong hipk = 0;
        for (ulong i = 0; i < pks.length; ++i)
        {
            if (nonPeriodic[pks[i]] > nonPeriodic[pks[hipk]])
                hipk = i;
        }

        // list of center peaks to test
        ulong[] testpk;
        testpk ~= hipk;
        if (hipk > 1)
            testpk ~= hipk - 1;
        if (hipk + 1 < pks.length)
            testpk ~= hipk + 1;

        //
        bool match = false;
        float[] periodic;
        for (ulong testidx = 0; !match && testidx < testpk.length; ++testidx)
        {
            periodic = extractPeriodicRange(nonPeriodic, samplerate, currentNote, cast(long)testpk[testidx]);

            float requiredPercent = 95;
            float actualPercent = 100.0f * periodic.length / nonPeriodic.length;
            stderr.writefln!"Attempt %d %f%%"(testidx + 1, actualPercent);
            match = actualPercent > requiredPercent;
        }

        if (!match)
            throw new Exception("Could not determine the periodic signal");

        seqs[iseq] = periodic;
    }

    /// crossfade samples
    for (ulong iseq = 0; iseq < seqs.length; ++iseq)
    {
        long pos = seqs[iseq].ptr - snd.channels[0].ptr;
        ulong len = seqs[iseq].length;
        enforce(pos >= 0 && pos + len <= snd.channels[0].length);
        seqs[iseq] = crossFadeSection(snd.channels[0], pos, pos + len);
    }

    /// save samples
    for (ulong iseq = 0; iseq < seqs.length; ++iseq)
    {
        string outputfile = outputPath ~ "/" ~ inputPath.baseName.stripExtension ~ "." ~ format!"%02d"(iseq) ~ ".wav";
        //stderr.writeln(outputfile);

        SF_INFO info = {};
        info.frames = seqs[iseq].length;
        info.channels = 1;
        ///info.format = SF_FORMAT_WAV|SF_FORMAT_FLOAT;
        info.format = SF_FORMAT_WAV|SF_FORMAT_PCM_16;
        info.samplerate = cast(int)samplerate;
        SNDFILE *sf = sf_open(outputfile.toStringz, SFM_WRITE, &info);
        if (!sf)
            throw new Exception("Cannot open wave file for writing");

        sf_writef_float(sf, seqs[iseq].ptr, seqs[iseq].length);

        scope(exit) sf_close(sf);
    }
}

MidiRange getMidiKeysOfFile(string inputPath)
{
    string name = inputPath.baseName.stripExtension;

    auto re = regex("(UK|LK|PK)_.*(4|8|16)");
    auto ma = name.match(re);

    if (!ma)
        throw new Exception("Cannot extract the midi key");

    MidiRange mr;
    if (ma.captures[1] == "UK")
    {
        mr.basekey = 53;
        mr.numkeys = 37;
    }
    else if (ma.captures[1] == "LK")
    {
        mr.basekey = 41;
        mr.numkeys = 37;
    }
    else if (ma.captures[1] == "PK")
    {
        mr.basekey = 36;
        mr.numkeys = 13;
    }
    else
        throw new Exception("Cannot extract the midi key");

    if (ma.captures[2] == "4")
        mr.basekey += 12;
    else if (ma.captures[2] == "8")
        mr.basekey += 0;
    else if (ma.captures[2] == "16")
        mr.basekey -= 12;
    else
        throw new Exception("Cannot extract the midi key");

    return mr;
}

float[] crossFadeSection(float[] wavedata, size_t start, size_t end, ulong xfsize = 16)
{
    if (start < xfsize || end + xfsize > wavedata.length)
        throw new Exception("Not enough samples");

    size_t count = end - start;
    float[] xfade = wavedata[start..end].dup;

    float[] rsrc = wavedata[start-xfsize..start];
    float[] lsrc = wavedata[end..end+xfsize];
    float[] ldst = xfade[0..xfsize];
    float[] rdst = xfade[count-xfsize..count];

    for (ulong i = 0; i < xfsize; ++i) {
        float mu = cast(float)(i + 1)/(xfsize + 1);
        //ldst[i] = mu * ldst[i] + (1-mu) * lsrc[i];
        rdst[i] = mu * rsrc[i] + (1-mu) * rdst[i];
    }

    return xfade;
}

float[] extractPeriodicRange(float[] wave, float samplerate, uint midiNote, long startPk = -1L)
{
    float currentFreq = 440.0f * pow(2.0f, (midiNote - 69.0f) / 12.0f);

    ulong[] pks = findPeaks(wave);
    if (pks.empty)
        throw new Exception("Wave does not have any peaks");

    // find highest peak
    ulong hipk = 0;
    if (startPk != -1L)
        hipk = cast(ulong)startPk;
    else
    {
        for (ulong i = 0; i < pks.length; ++i)
        {
            if (wave[pks[i]] > wave[pks[hipk]])
                hipk = i;
        }
        //stderr.writefln!"hipeak=%f"(pks[hipk] / samplerate);
    }

    // find periods left
    ulong pk1 = hipk;
    for (ulong i = hipk; i-- > 0; )
    {
        float thisPeriod = (pks[pk1] - pks[i]) / samplerate;
        float thisFreq = 1.0f / thisPeriod;
        float thisMidiPitch = 12.0f * log2(thisFreq / 440.0f) + 69.0f;
        uint thisMidiNote = cast(uint)lround(thisMidiPitch);

        bool validPeriod =
            (thisMidiNote == midiNote
             || cast(int)thisMidiNote == cast(int)midiNote - 12
             || cast(int)thisMidiNote == cast(int)midiNote - 24
             // || cast(int)thisMidiNote == cast(int)midiNote - 36
             // || cast(int)thisMidiNote == cast(int)midiNote - 48
             // || cast(int)thisMidiNote == cast(int)midiNote - 60
             // || cast(int)thisMidiNote == cast(int)midiNote - 72
             // || cast(int)thisMidiNote == cast(int)midiNote - 96
            );

        //stderr.writefln!"[L] note=%f freq=%f period=%f timepos=%f valid=%d"(thisMidiPitch, thisFreq, thisPeriod, pks[i] / samplerate, validPeriod);

        if (validPeriod)
            pk1 = i;
    }
    // find periods right
    ulong pk2 = hipk;
    for (ulong i = hipk + 1; i < pks.length; ++i)
    {
        float thisPeriod = (pks[i] - pks[pk2]) / samplerate;
        float thisFreq = 1.0f / thisPeriod;
        float thisMidiPitch = 12.0f * log2(thisFreq / 440.0f) + 69.0f;
        uint thisMidiNote = cast(uint)lround(thisMidiPitch);

        bool validPeriod =
            (thisMidiNote == midiNote
             || cast(int)thisMidiNote == cast(int)midiNote - 12
             || cast(int)thisMidiNote == cast(int)midiNote - 24
             // || cast(int)thisMidiNote == cast(int)midiNote - 36
             // || cast(int)thisMidiNote == cast(int)midiNote - 48
             // || cast(int)thisMidiNote == cast(int)midiNote - 60
             // || cast(int)thisMidiNote == cast(int)midiNote - 72
             // || cast(int)thisMidiNote == cast(int)midiNote - 84
             // || cast(int)thisMidiNote == cast(int)midiNote - 96
            );

        //stderr.writefln!"[R] note=%f freq=%f period=%f timepos=%f valid=%d"(thisMidiPitch, thisFreq, thisPeriod, pks[i] / samplerate, validPeriod);

        if (validPeriod)
            pk2 = i;
    }

    return wave[pks[pk1]..pks[pk2]];
}

ulong[] findPeaks(float[] wave, float maxRelErrToHighest = 0.2f, ulong neighboors = 10)
{
    ulong[] pks = [];
    pks.reserve(256);

    ulong count = wave.length;
    if (count == 0)
        return pks;

    /// find highest peak
    ulong highest = 0;
    for (ulong i = 1; i < count; ++i) {
        if (abs(wave[i]) > abs(wave[highest]))
            highest = i;
    }
    int highsign = signbit(wave[highest]);

    /// find local peaks within a set neighborhood
    /// having the same polarity as the highest
    for (ulong i = neighboors; i + neighboors < count; ++i) {
        bool ispeak = true;

        if (signbit(wave[i]) != highsign)
            ispeak = false;

        for (ulong j = 1; j < neighboors && ispeak; ++j) {
            // XXX equal case not supported. probably not important
            float y = wave[i];
            float y1 = wave[i+j];
            float y2 = wave[i-j];
            ispeak =
                (signbit(y) != signbit(y1) || abs(y) > abs(y1)) &&
                (signbit(y) != signbit(y2) || abs(y) > abs(y2));
        }

        // keep peaks only within some % of highest
        if (ispeak) {
            float rel = abs(wave[i] - wave[highest]) / abs(wave[highest]);
            //stderr.writefln!"PEAK cur=%f high=%f rel=%f"(wave[i], wave[highest], rel);
            if (rel < maxRelErrToHighest)
                pks ~= i;
        }
    }

    return pks;
}

float[][] reduceToCenterSamples(float[][] seqs, float srate, float duration)
{
    ulong nseq = seqs.length;
    float[][] newSeqs = new float[][nseq];

    ulong durframes = cast(ulong)ceil(duration * srate);

    for (ulong i = 0; i < nseq; ++i)
    {
        if (seqs[i].length < durframes)
            throw new Exception("Not enough samples");
        ulong offset = (seqs[i].length - durframes) / 2;
        newSeqs[i] = seqs[i][offset..offset+durframes];
    }

    return newSeqs;
}

float[][] extractNonZeroSequences(float[] data, float srate, float tmin, float tol = 1e-3f)
{
    float[][] seqs;
    ulong total = data.length;

    ulong numthr = cast(ulong)ceil(tmin * srate);

    seqs.reserve(256);

    ulong i = 0;
    while (i < total && abs(data[i]) < tol)
        i += 1;

    ulong start = i;
    while (i < total)
    {
        ulong cntz = 0;
        while (cntz < total - i && abs(data[i + cntz]) < tol)
            ++cntz;
        if (cntz == 0)
            ++i;
        else
        {
            if (cntz < numthr)
                i += cntz;
            else
            {
                seqs ~= data[start..i];
                i += cntz;
                start = i;
            }
        }
    }

    return seqs;
}

Sound *readSoundFile(string filename)
{
    SF_INFO info;
    SNDFILE *sf = sf_open(filename.toStringz, SFM_READ, &info);
    if (!sf)
        throw new Exception("Cannot open the sound file");
    scope(exit) sf_close(sf);

    Sound *snd = new Sound;
    ulong numFrames = cast(ulong)info.frames;
    ulong numChannels = cast(ulong)info.channels;
    snd.samplerate = info.samplerate;
    snd.frames = numFrames;
    snd.channels = new float[][numChannels];

    for (ulong c = 0; c < numChannels; ++c)
        snd.channels[c] = new float[snd.frames];

    float[] interleaved = new float[numFrames * numChannels];
    sf_readf_float(sf, interleaved.ptr, numFrames);
    for (ulong c = 0; c < numChannels; ++c)
    {
        for (ulong i = 0; i < numFrames; ++i)
            snd.channels[c][i] = interleaved[i * numChannels + c];
    }

    return snd;
}
