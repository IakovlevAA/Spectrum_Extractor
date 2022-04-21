import xml.etree.ElementTree as ET
import numpy as np
import csv
from scipy.signal import find_peaks

def write_csv_first(fieldnames, data):
    with open('Results.csv', 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()

        for f in range(len(data[0])):
            dict = {}
            dict[fieldnames[0]] = data[0][f]
            dict[fieldnames[3]] = data[3][f]
            dict[fieldnames[4]] = data[4][f]
            dict[fieldnames[5]] = data[5][f]

            for f_q in range(len(data[1])):
                dict[fieldnames[1]] = data[1][f_q]
                dict[fieldnames[2]] = data[2][f_q]

                for n in range(0, len(data[6][f_q])):
                    dict[fieldnames[n + 6]] = data[6][f_q][n]
                writer.writerow(dict)


def write_csv_second(fieldnames, data):
    with open('Results2.csv', 'w') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
        writer.writeheader()

        for f in range(len(data[0])):
            dict = {}
            dict[fieldnames[0]] = data[0][f]
            for f_q in range(len(data[1])):
                dict[fieldnames[1]] = data[1][f_q]
                dict[fieldnames[2]] = data[2][f_q]
                dict[fieldnames[9]] = data[4][f_q]
                dict[fieldnames[10]] = data[5][f_q]
                dict[fieldnames[11]] = data[6][f_q]
                for n in range(0, len(data[3][f_q])):
                    dict[fieldnames[n + 3]] = data[3][f_q][n]
                writer.writerow(dict)


def samples_for_first(FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate):
    is_channel = []
    N = 1024 * 1024
    sample_table = []
    for i in range(len(FilePath)):
        # samples = np.fromfile(FilePath[i], dtype=np.short, count=2 * N)
        samples = np.fromfile('04.06.2015 18_57_07_1212,98MHz 20722.2KHz.pcm', dtype=np.short, count=2 * N)
        samples = samples[0::2] + 1j * samples[1::2]
        samples = np.fft.fftshift(np.fft.fft(samples))


        for j in range(len(freq_center_signal)):
            start_pos = int(N * (freq_center_signal[j] - bandwidth_signal[j] / 2) / sample_rate[i] + N / 2)
            end_pos = int(N * (freq_center_signal[j] + bandwidth_signal[j] / 2) / sample_rate[i] + N / 2)
            s = samples[start_pos:end_pos]

            # значения для таблицы

            is_channel.append(1)

            s = 10 * np.log10(np.abs(s))
            sample_table.append(np.resize(s, 64))
    fields = ['FilePath', 'freq_center_signal', 'bandwidth_signal', 'freq_center', 'sample_rate', 'is_channel']
    for i in range(0, len(sample_table[0])):
        fields.append('samples' + str(i + 1))

    data = [FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate, is_channel, sample_table]
    return data, fields


def samples_for_second(FilePath, freq_center_signal, bandwidth_signal, freq_center_record, sample_rate,modulation):
    fft_size = 1024 * 1024
    pwr_table = []
    env_table = []
    freq_table = []
    for i in range(len(FilePath)):
        # samples = np.fromfile(FilePath[i], dtype=np.short, count=2 * fft_size)
        samples = np.fromfile('04.06.2015 18_57_07_1212,98MHz 20722.2KHz.pcm', dtype=np.short, count=2 * fft_size)
        samples = samples / max(samples)
        samples = samples[0::2] + 1j * samples[1::2]

        for j in range(len(freq_center_signal)):

            freq_offs_rel = (freq_center_signal[j] - freq_center_record[i]) / sample_rate[i]
            phases = np.linspace(0, len(samples) - 1, len(samples))
            samples = samples * np.exp(-1j * 2 * np.pi * freq_offs_rel * phases)

            pwrs = [2, 4, 6, 8, 12, 16]
            samples = samples / max(samples)
            # pwr
            p = []
            for pwr_ctr in range(len(pwrs)):
                pwr = pwrs[pwr_ctr]
                samples_freq = np.fft.fftshift(np.fft.fft(samples[0:fft_size] ** pwr))
                samples_freq_abs = np.abs(samples_freq)
                samples_freq_abs = 10.0 * np.log10(samples_freq_abs + 1.0E-10)
                peaks = find_peaks(samples_freq_abs)
                p.append(len(peaks[0]))
            pwr_table.append(p)
            # env
            samples_freq = np.fft.fftshift(np.fft.fft(np.abs(samples[0:fft_size])))
            samples_freq_abs = np.abs(samples_freq)

            peaks = find_peaks(samples_freq_abs)
            env_table.append(len(peaks[0]))
            # freq
            samples_dif_angle = np.angle(samples[1:fft_size] * np.conj(samples[0:fft_size - 1]))
            samples_freq = np.fft.fftshift(np.fft.fft(np.abs(samples_dif_angle)))
            samples_freq_abs = np.abs(samples_freq)

            peaks = find_peaks(samples_freq_abs)
            freq_table.append(len(peaks[0]))
    fields = ['FilePath', 'freq_center_signal', 'bandwidth_signal', 'pwr2','pwr4','pwr6','pwr8','pwr12','pwr16', 'env', 'freq','modulation']
    data = [FilePath, freq_center_signal, bandwidth_signal, pwr_table, env_table, freq_table, modulation]
    return data, fields

def parse(xmlFile):
    FilePath = []
    freq_center_signal = []
    bandwidth_signal = []
    freq_center = []
    sample_rate = []
    modulation=[]
    root_node = ET.parse(xmlFile).getroot()
    for tag in root_node.findall('Test/SignalsBunch/SignalStream/FilePath'):
        FilePath.append(tag.text)

    for tag in root_node.findall('Test/SignalsInSelection/Signal/Location'):
        freq_center_signal.append(float(tag.attrib['freq_center']))
        bandwidth_signal.append(float(tag.attrib['bandwidth']))

    for tag in root_node.findall('Test/SignalsBunch/SignalStream/StreamParams'):
        freq_center.append(float(tag.attrib['freq_center']))
        sample_rate.append(float(tag.attrib['sample_rate']))

    for tag in root_node.findall('Test/SignalsInSelection/Signal/SignalType'):
        modulation.append(float(tag.attrib['modulation']))


    return FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate,modulation


FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate,modulation = parse('TestLocalizerTest.xml')
data, fields = samples_for_first(FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate)
write_csv_first(fields, data)

data, fields = samples_for_second(FilePath, freq_center_signal, bandwidth_signal, freq_center, sample_rate,modulation)
write_csv_second(fields, data)