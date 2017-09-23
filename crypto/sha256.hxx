#include "detail/sha256.hxx"

namespace crypto {

template <class word_t, class byte_t> void sha256(word_t h[8], const byte_t *data, size_t n)
{
	h[0] = 0x6a09e667;	h[1] = 0xbb67ae85;	h[2] = 0x3c6ef372;	h[3] = 0xa54ff53a;
	h[4] = 0x510e527f;	h[5] = 0x9b05688c;	h[6] = 0x1f83d9ab;	h[7] = 0x5be0cd19;
	detail::chunks_stream<word_t, byte_t> stream(data, n);
    for (word_t w[64]; stream.next(w); ) {
    	for (int i = 16; i < 64; ++i) {
    		w[i] = w[i - 16] + w[i - 7] 
    			+ (detail::rots<7, 18>(w[i - 15]) ^ (w[i - 15] >> 3))
    			+ (detail::rots<17, 19>(w[i - 2]) ^ (w[i - 2] >> 10));
    	}
		word_t t[8] = {h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]};
		detail::sha256_helper<0>::round(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], w);
		for (int i = 0; i < 8; ++i)
			h[i] += t[i];
	}
}

template <class word_t, class Data> void sha256(word_t h[8], const Data &data)
{
	sha256(h, data.data(), data.size());
}

} // namespace crypto