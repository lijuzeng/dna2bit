//
// Created by ljz on 23-12-29.
//

#ifndef DNA_TO_BIT_DISTANCE_HPP
#define DNA_TO_BIT_DISTANCE_HPP
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <omp.h>
using namespace std;

void LoadBit(const char *file_name, vector<uint64_t> &bit) {
    FILE *f = fopen(file_name, "rb");
    if (f == nullptr) {cerr << "File open fail" << endl;    exit(0);}
    uint64_t tmp = 0;
    while (fread(&tmp, 8, 1, f) == 1) {bit.push_back(tmp);}
    fclose(f);
}

float distance_metric(const vector<uint64_t> &bit1, const vector<uint64_t> &bit2, const size_t kmer, const size_t bit) {
    size_t len1, len2;
    float cnt = 0;
    len1 = bit1.size();     len2 = bit2.size();
    if (len1 != len2) {cerr << "two bits of differen lenght" << endl;   exit(0);}
    for (size_t i = 0; i < len1; ++i) cnt += (float)__builtin_popcountll(bit1[i] ^ bit2[i]);
    cnt = -1 / (float)kmer * log(1 - cnt / (float)bit);

    return cnt;
}

void triangular_distance_matrix(vector<string> &file_names_list, const string &output_filename, const size_t threads, const size_t kmer, const size_t bit) {
    FILE *out = fopen(output_filename.c_str(), "w");
    vector<vector<uint64_t>> bit_list;
    for (const auto &i : file_names_list) {
        vector<uint64_t> tmp;
        LoadBit(i.c_str(), tmp);
        bit_list.push_back(tmp);
    }
    size_t length = bit_list.size();
    fputs("sample", out);
    for (const auto &i : file_names_list) fprintf(out, ",%s", i.c_str());
    fputc('\n', out);
#pragma omp parallel for schedule(static, 1) num_threads(threads ? threads : 1)
    for (size_t i = 0; i < length; ++i) {
        FILE *lineout = fopen((to_string(i) + ".txt").c_str(), "w");
        fputs(file_names_list[i].c_str(),lineout);
        for (size_t j = 0; j < i; ++j) fputs(",-", lineout);
        fputs(",0", lineout);
        for (size_t j = i + 1;j < length; ++j) fprintf(lineout, ",%f", distance_metric(bit_list[i],bit_list[j], kmer, bit));
        fputc('\n', lineout);
        fclose(lineout);
    }
    for(size_t i = 0; i < length; ++i){
        string fin = to_string(i) + ".txt";
        FILE *linein = fopen(fin.c_str(), "r");
        char buftmp[1024];
        while (fgets(buftmp, 1024, linein) != nullptr) {fputs(buftmp,out);}
        fclose(linein);
        remove(fin.c_str());
    }
    fclose(out);
}

void whole_distance_matrix(vector<string> &file_names_list, const string &output_filename, const size_t threads, const size_t kmer, const size_t bit) {
    FILE *out = fopen(output_filename.c_str(), "w");
    vector<vector<uint64_t>> bit_list;
    for (const auto &i : file_names_list) {
        vector<uint64_t> tmp;
        LoadBit(i.c_str(), tmp);
        bit_list.push_back(tmp);
    }
    size_t length = bit_list.size();
    fputs("sample", out);
    for (const auto &i : file_names_list) fprintf(out, ",%s", i.c_str());
    fputc('\n', out);
#pragma omp parallel for schedule(static, 1) num_threads(threads ? threads : 1)
    for (size_t i = 0; i < length; ++i) {
        FILE *lineout = fopen((to_string(i) + ".txt").c_str(), "w");
        fputs(file_names_list[i].c_str(),lineout);
        for (size_t j = 0; j < length; ++j) fprintf(lineout, ",%f", distance_metric(bit_list[i],bit_list[j], kmer, bit));
        fputc('\n', lineout);
        fclose(lineout);
    }
    for(size_t i = 0; i < length; ++i){
        string fin = to_string(i) + ".txt";
        FILE *linein = fopen(fin.c_str(), "r");
        char buftmp[1024];
        while (fgets(buftmp, 1024, linein) != nullptr) {fputs(buftmp,out);}
        fclose(linein);
        remove(fin.c_str());
    }
    fclose(out);
}

#endif //DNA_TO_BIT_DISTANCE_HPP
