//
// Created by tyx and ljz on 24-1-8.
//

#ifndef DNA_TO_BIT_ROLLINGHASH_H
#define DNA_TO_BIT_ROLLINGHASH_H
#include <cstdint>
#include <iostream>

const uint8_t doff = 0x07;
static const uint64_t hA = 0x3c8bfbb395c60474,hC = 0x3193c18562a02b4c,hG = 0x20323ed082572324,hT = 0x295549f54be24456;
static const uint64_t Tab[256] = {0,hT,0,hG,hA,0,0,hC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,hA,0,hC,0,0,0,hG,0,0,0,0,0,0,0,0,0,0,0,0,hT,0,0,0,0,0,0,0,0,0,0,0,0,hA,0,hC,0,0,0,hG,0,0,0,0,0,0,0,0,0,0,0,0,hT,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

inline uint64_t r33(const uint64_t vv) {
    uint64_t v=(vv << 1) | (vv >> 63);
    uint64_t x = (v ^ (v >> 33)) & 1;
    return v ^ (x | (x << 33));
}

inline uint64_t r3263(const uint64_t vv) {
    uint64_t v=(vv >> 1) | (vv << 63);
    uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
}

inline uint64_t rol(unsigned b, const uint64_t v, unsigned s)  {
    return ((v << s%b) | (v >> (b - s%b)));
}

inline uint64_t opchar(const unsigned char ch, const size_t k) {
    uint64_t tmp1=(rol(31,Tab[ch]>>33,k)&0x7FFFFFFF)<<33;
    uint64_t tmp2=rol(33,Tab[ch]&0x1FFFFFFFF,k)& 0x1FFFFFFFF;
    return tmp1 | tmp2;
}
#endif //DNA_TO_BIT_ROLLINGHASH_H
