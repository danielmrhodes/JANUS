#ifndef DATA_FORMAT_h
#define DATA_FORMAT_h 1

struct Header {
  
  int evtNum, nBdata, nSdata;

  size_t bytes() {
    return 3*sizeof(int);
  }
  
}__attribute__((__packed__));

struct SegaData {

  int det, seg;
  double en, x, y, z;
  bool fep, pfep;
  
};

struct Bambino2Data {

  int det, ring, sector;
  double en, x, y, z;
  bool proj, rec;

  bool IsRing() const {return (bool)ring;}

  bool operator<(const Bambino2Data& rhs) const { return en > rhs.en; }
  
};

struct JANUSData {

  Bambino2Data bData[10]; //100
  SegaData sData[50]; //100

}__attribute__((__packed__));

#endif
