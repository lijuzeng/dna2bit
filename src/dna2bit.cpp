//
// Created by ljz on 23-12-05.
// VERSION 0.1.8
//
#include <iostream>
#include "hashmain.hpp"
#include "command.hpp"
#include "distance.hpp"
#include <omp.h>
using namespace std;

int main(int argc, char *argv[]) {
    char code[256];
    memset(code, 'N', 256);
    code['A'] = 'T';
    code['T'] = 'A';
    code['G'] = 'C';
    code['C'] = 'G';
    if (argc < 2 || !strcmp(argv[1], "--help")) { MainUsage(argv);}
    if (!strcmp(argv[1], "sketch")) {
        SketchArgs sket;
        sket.GetArgs(argc, argv);
        if (sket.list_files_name.empty()) {
            cerr << argv[optind + 1] << "\t encoding ..." << endl;
            string out_file = (string)argv[optind + 1] + ".k." + to_string(sket.kmer_len) + ".l." + to_string(sket.bit_len) + ".bit";
            File2Bit(argv[optind + 1], out_file.c_str(), code, sket.bit_len, sket.kmer_len, sket.hash_type);
            cerr << argv[optind + 1] << "\t encoded !!!" << endl;
        } else {
            size_t files_size = sket.file_names_list.size();
            cerr << "files of " << sket.list_files_name << " encoding..." << endl;
            cerr << "\tusing " << (sket.nthreads ? sket.nthreads : omp_get_max_threads()) << " threads" << endl;
#pragma omp parallel for schedule(static, 1) num_threads(sket.nthreads ? sket.nthreads : omp_get_max_threads())
            for (size_t i = 0; i < files_size; ++i) {
                File2Bit(sket.file_names_list[i].c_str(), (sket.file_names_list[i] + ".k." + to_string(sket.kmer_len) + ".l." +
                        to_string(sket.bit_len) + ".bit").c_str(), code, sket.bit_len, sket.kmer_len, sket.hash_type);
            }
            cerr << "files of " << sket.list_files_name << " encoded!!!" << endl;
        }
    } else if (!strcmp(argv[1], "dist")) {
        DistArgs dist;
        dist.GetArgs(argc, argv);
        if (dist.list_files_name.empty()) {
            vector<uint64_t> bit0, bit1;
            LoadBit(argv[optind + 1], bit0);    LoadBit(argv[optind + 2], bit1);
            cout << "distance metric : " << distance_metric(bit0, bit1, dist.kmer_len, dist.bit_len) << endl;
        } else {
            if (dist.triangular) {
                cerr << "triangular distance matrix calculating ..." << endl;
                cerr << "\tusing " << (dist.nthreads ? dist.nthreads : 1) << " threads" << endl;
                triangular_distance_matrix(dist.file_names_list, dist.output_filename, dist.nthreads, dist.kmer_len, dist.bit_len);
                cerr << "triangular distance matrix calculated !!!" << endl;
            } else {
                cerr << "whole distance matrix calculating ..." << endl;
                cerr << "\tusing " << (dist.nthreads ? dist.nthreads : 1) << " threads" << endl;
                whole_distance_matrix(dist.file_names_list, dist.output_filename, dist.nthreads, dist.kmer_len, dist.bit_len);
                cerr << "whole distance matrix calculated !!!" << endl;
            }
        }
    } else {
        cerr << "Unrecognized command : " << argv[1] << endl;
        MainUsage(argv);
    }
    return 0;
}