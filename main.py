import xml.etree.ElementTree as ET
import numpy as np
import csv
from scipy.signal import find_peaks


def write_csv_first(fieldnames, data):
    with open('Results.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()

        for f in range(len(data)):
            writer.writerow(data[f])


def write_csv_second(fieldnames, data):
    with open('Results2.csv', 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()

        for f in range(len(data)):
            writer.writerow(data[f])


def samples_for_first(data):
    N = 1024 * 1024
    dicts = []
    for i in range(len(data)):
        dict = {}
        dict = {'FilePath': data[i]['FilePath'], 'freq_center_signal': data[i]['freq_center_signal'],
                'bandwidth_signal': data[i]['bandwidth_signal'],
                'freq_center': data[i]['freq_center'], 'sample_rate': data[i]['sample_rate'], 'is_channel': 1}

        # samples = np.fromfile(FilePath[i], dtype=np.short, count=2 * N)
        samples = np.fromfile('04.06.2015 18_57_07_1212,98MHz 20722.2KHz.pcm', dtype=np.short, count=2 * N)
        samples = samples[0::2] + 1j * samples[1::2]
        samples = np.fft.fftshift(np.fft.fft(samples))

        start_pos = int(
            N * (data[i]['freq_center_signal'] - data[i]['bandwidth_signal'] / 2) / data[i]['sample_rate'] + N / 2)
        end_pos = int(
            N * (data[i]['freq_center_signal'] + data[i]['bandwidth_signal'] / 2) / data[i]['sample_rate'] + N / 2)

        s = samples[start_pos:end_pos]
        s = 10 * np.log10(np.abs(s))
        s = np.resize(s, 64)

        for j in range(0, len(s)):
            dict['sample' + str(j + 1)] = s[j]
        dicts.append(dict)

    fields = ['FilePath', 'freq_center_signal', 'bandwidth_signal', 'freq_center', 'sample_rate', 'is_channel']
    for i in range(0, 64):
        fields.append('sample' + str(i + 1))

    return dicts, fields


def samples_for_second(data):
    fft_size = 1024 * 1024
    dicts = []
    pwrs = [2, 4, 6, 8, 12, 16]
    for i in range(len(data)):
        dict = {}
        dict = {'FilePath': data[i]['FilePath'], 'freq_center_signal': data[i]['freq_center_signal'],
                'bandwidth_signal': data[i]['bandwidth_signal']}
        # samples = np.fromfile(FilePath[i], dtype=np.short, count=2 * fft_size)
        samples = np.fromfile('04.06.2015 18_57_07_1212,98MHz 20722.2KHz.pcm', dtype=np.short, count=2 * fft_size)
        samples = samples / max(samples)
        samples = samples[0::2] + 1j * samples[1::2]

        freq_offs_rel = (data[i]['freq_center_signal'] - data[i]['freq_center']) / data[i]['sample_rate']
        phases = np.linspace(0, len(samples) - 1, len(samples))
        samples = samples * np.exp(-1j * 2 * np.pi * freq_offs_rel * phases)

        samples = samples / max(samples)
        # pwr
        for pwr_ctr in range(len(pwrs)):
            pwr = pwrs[pwr_ctr]
            samples_freq = np.fft.fftshift(np.fft.fft(samples[0:fft_size] ** pwr))
            samples_freq_abs = np.abs(samples_freq)
            samples_freq_abs = 10.0 * np.log10(samples_freq_abs + 1.0E-10)
            peaks = find_peaks(samples_freq_abs)
            dict['pwr' + str(pwr)] = len(peaks[0])

        # env
        samples_freq = np.fft.fftshift(np.fft.fft(np.abs(samples[0:fft_size])))
        samples_freq_abs = np.abs(samples_freq)

        peaks = find_peaks(samples_freq_abs)
        dict['env'] = len(peaks[0])
        # freq
        samples_dif_angle = np.angle(samples[1:fft_size] * np.conj(samples[0:fft_size - 1]))
        samples_freq = np.fft.fftshift(np.fft.fft(np.abs(samples_dif_angle)))
        samples_freq_abs = np.abs(samples_freq)

        peaks = find_peaks(samples_freq_abs)
        dict['freq'] = len(peaks[0])

        dict['modulation'] = data[i]['modulation']
        dicts.append(dict)
    fields = ['FilePath', 'freq_center_signal', 'bandwidth_signal', 'pwr2', 'pwr4', 'pwr6', 'pwr8', 'pwr12', 'pwr16',
              'env', 'freq', 'modulation']

    return dicts, fields


def parse(xmlFile):
    root_node = ET.parse(xmlFile).getroot()
    dicts = []
    for element in root_node.iterfind('Test'):
        for el1 in element.iterfind('SignalsInSelection/Signal'):
            data = {}
            data['FilePath'] = element.find('SignalsBunch/SignalStream/FilePath').text
            data['freq_center'] = float(element.find('SignalsBunch/SignalStream/StreamParams').attrib['freq_center'])
            data['sample_rate'] = float(element.find('SignalsBunch/SignalStream/StreamParams').attrib['sample_rate'])
            data['freq_center_signal'] = float(el1.find('Location').attrib['freq_center'])
            data['bandwidth_signal'] = float(el1.find('Location').attrib['bandwidth'])
            data['modulation'] = float(el1.find('SignalType').attrib['modulation'])
            dicts.append(data)
    return dicts


data = parse('TestLocalizerTest.xml')
data1, fields = samples_for_first(data)
write_csv_first(fields, data1)

data2, fields = samples_for_second(data)
write_csv_second(fields, data2)
