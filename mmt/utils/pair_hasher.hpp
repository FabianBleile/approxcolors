#ifndef PAIR_HASHER_H_
#define PAIR_HASHER_H_

struct UInt32PairHash {
  std::size_t operator()(const std::pair<uint32_t, uint32_t> &p) const  {
    assert(sizeof(std::size_t)>=8);  //Ensure that std::size_t, the type of the hash, is large enough
    //Shift first integer over to make room for the second integer. The two are
    //then packed side by side.
    return (((uint64_t)p.first)<<32) | ((uint64_t)p.second);
  }
};

#endif
