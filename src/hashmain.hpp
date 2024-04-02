//
// Created by ljz on 24-1-5.
//

#ifndef DNA_TO_BIT_HASHMAIN_HPP
#define DNA_TO_BIT_HASHMAIN_HPP
#include "wyhash.h"
#include "rollinghash.h"
#include "MurmurHash3.h"
#include <unistd.h>
#include "zlib.h"
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include "kseq.h"
using namespace std;
KSEQ_INIT(gzFile, gzread)

/* read 序列编码成 distance vector */
void read2dis_wy(const char *read, const size_t sliding_len, vector<long> &dis) {
    size_t length = strlen(read), vec_size = dis.size();
    for (size_t i = 0; i < length - sliding_len; ++i) {
        uint64_t hash = wyhash(read + i, sliding_len, 0, _wyp);
        dis[(hash & 0x7FFFFFFFFFFFFFFF) % (vec_size)] += 2 * (long)(hash >> 63) - 1;
    }
}
void read2dis_mu(const char *read, const size_t sliding_len, vector<long> &dis) {
    size_t length = strlen(read), vec_size = dis.size();
    char data[16];
    uint64_t hash;
    for (size_t i = 0; i < length - sliding_len; ++i) {
        MurmurHash3_x64_128(read, sliding_len, 0, data);
        hash = *((uint64_t *) data);
        dis[(hash & 0x7FFFFFFFFFFFFFFF) % (vec_size)] += 2 * (long)(hash >> 63) - 1;
    }
}
/* read 与其互补链的 distance vector */
void reads2dis_wy(const char *code, const char *read, const size_t sliding_len, vector<long> &dis) {
    char *com;
    size_t length = strlen(read);
    com = (char *)malloc(length * sizeof(char));
    for (size_t i = 0; i < length - 1; ++i) com[i] = code[(int)read[length - 2 - i]];
    com[length - 1] = '\n';
    read2dis_wy(read, sliding_len, dis);
    read2dis_wy(com, sliding_len, dis);
    free(com);
}
void reads2dis_mu(const char *code, const char *read, const size_t sliding_len, vector<long> &dis) {
    char *com;
    size_t length = strlen(read);
    com = (char *)malloc(length * sizeof(char));
    for (size_t i = 0; i < length - 1; ++i) com[i] = code[(int)read[length - 2 - i]];
    com[length - 1] = '\n';
    read2dis_mu(read, sliding_len, dis);
    read2dis_mu(com, sliding_len, dis);
    free(com);
}
/* 将 distance vector 变为 bit vector */
inline void dis2bit(vector<long> &bit) {
    size_t len = bit.size();
    for (size_t i = 0; i < len; ++i) {bit[i] = (bit[i] >> 63) + 1;}
}
/* 对 基因组 文件的处理 */
void HandleFaFile_wy(gzFile *in, const char *code, size_t sli, vector<long> &bit) {
    kseq_t *seq = kseq_init(*in);
    while (kseq_read(seq) >= 0) reads2dis_wy(code, seq->seq.s, sli, bit);
    kseq_destroy(seq);
}
void HandleFaFile_ro(gzFile *in, size_t sli, vector<long> &bit) {
    size_t len = bit.size();
    kseq_t *seq = kseq_init(*in);
    while (kseq_read(seq) >= 0) {
        uint64_t f = 0,r = 0;
        string read = seq->seq.s;
        for (int i = sli - 1; i >= 0; --i) {
            f = r33(f) ^ Tab[(int)read[sli - 1 - i]];
            r = r33(r) ^ Tab[read[i] & doff];
        }
        bit[(f & 0x7FFFFFFFFFFFFFFF) % (len)] += (2 * int(f >> 63) - 1);
        bit[(r & 0x7FFFFFFFFFFFFFFF) % (len)] += (2 * int(r >> 63) - 1);
        for (size_t i = 0; i < read.length() - sli; ++i){
            f = r33(f) ^ Tab[(int)read[i + sli]] ^ opchar(read[i],sli);
            r = r3263(r ^ opchar(read[i + sli] & doff,sli)^Tab[read[i] & doff]);
            bit[(f & 0x7FFFFFFFFFFFFFFF) % (len)] += (2 * int(f >> 63) - 1);
            bit[(r & 0x7FFFFFFFFFFFFFFF) % (len)] += (2 * int(r >> 63) - 1);
        }
    }
    kseq_destroy(seq);
}
void HandleFaFile_mu(gzFile *in, const char *code, size_t sli, vector<long> &bit) {
    kseq_t *seq = kseq_init(*in);
    while (kseq_read(seq) >= 0) reads2dis_mu(code, seq->seq.s, sli, bit);
    kseq_destroy(seq);
}
/* 写出 二进制 的 bit 文件 */
void BitOut(FILE *out, const vector<long> &bit) {
    size_t len = bit.size(), cnt = len / 64;
    for (size_t i = 0; i < cnt; ++i) {
        size_t loc = i * 64;
        uint64_t result = 0;
        for (size_t j = 0; j < 64; ++j) result = (result << 1) | bit[loc + j];
        fwrite(&result, 8, 1, out);
    }
    if (len % 64) {
        uint64_t result = 0;
        for (size_t i = cnt * 64; i < len; ++i) result = (result << 1) | bit[i];
        fwrite(&result, 8, 1, out);
    }

}

void File2Bit(const char *infile, const char *outfile, const char *code, size_t len, size_t sli, size_t hat) {
    vector<long> bit(len, 0);
    FILE *o = fopen(outfile, "wb");
    gzFile f = gzopen(infile, "r");
    if (f == nullptr) {
        printf("%s unable open\n", outfile);
        exit(0);
    }
    switch (hat) {
        case 0 : HandleFaFile_wy(&f, code, sli, bit);           break;
        case 1 : HandleFaFile_ro(&f, sli, bit);                 break;
        case 2 : HandleFaFile_mu(&f, code, sli, bit);           break;
        default : cerr << "Unrecognized hash type : " << hat << endl;   break;
    }
    dis2bit(bit);
    BitOut(o, bit);
    fclose(o);
    gzclose(f);
}
#endif //DNA_TO_BIT_HASHMAIN_HPP
