import xml.etree.ElementTree as ET
import numpy as np
import csv


def write_csv(fieldnames, data):
    with open('Results.csv', 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()
        dict = {}
        for n in range(0,len(fieldnames)):
            if 'samples' in fieldnames[n]:
                dict[fieldnames[n]]=data[len(data)-1][n-(len(data)-1)]
            else:
                dict[fieldnames[n]] = data[n]
        writer.writerow(dict)


def parse(xmlFile):
    root_node = ET.parse(xmlFile).getroot()
    for tag in root_node.findall('Test/SignalsBunch/SignalStream/FilePath'):
        FilePath = tag.text
    for tag in root_node.findall('Test/SignalsBunch/SignalStream/Selection'):
        freq_start = float(tag.attrib['freq_start'])
        freq_end = float(tag.attrib['freq_end'])
    for tag in root_node.findall('Test/SignalsBunch/SignalStream/StreamParams'):
        freq_center = float(tag.attrib['freq_center'])
        sample_rate = float(tag.attrib['sample_rate'])
    return FilePath, freq_start, freq_end, freq_center, sample_rate


if __name__ == '__main__':
    FilePath, freq_start, freq_end, freq_center, sample_rate = parse('TestLocalizerTest.xml')
    is_channel=1
    N = 1024 * 1024
    samples = np.fromfile('04.06.2015 18_57_07_1212,98MHz 20722.2KHz.pcm', dtype=np.short, count=2 * N)
    samples = samples[0::2] + 1j * samples[1::2]
    samples = np.fft.fftshift(np.fft.fft(samples))

    start_pos = int(N * (freq_start - freq_center) / sample_rate + N / 2)
    end_pos = int(N * (freq_end - freq_center) / sample_rate + N / 2)
    samples = samples[start_pos:end_pos]
    samples = np.resize(samples, 64)

    #модуль и децибелы
    abs_samples = np.abs(samples)
    decibels = np.log10(samples)

    fields = ['FilePath', 'freq_start', 'freq_end', 'freq_center', 'sample_rate','is_channel']
    for i in range(0, len(samples)):
        fields.append('samples' + str(i+1))
    data = [FilePath, freq_start, freq_end, freq_center, sample_rate,is_channel, samples]
    write_csv(fields, data)