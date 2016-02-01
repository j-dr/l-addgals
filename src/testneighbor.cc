int main(int argc, char* argv[])
{
  int id, ThisBCG;
  float mag, dens, ThisZ;
  vector <GalSED> galseds = ReadSED();
  mag = -12.315;
  dens = -0.001127588;
  ThisZ = 0.11;
  ThisBCG = 0;
  id = findCloseGalaxies2(galseds, mag, dens, ThisZ, ThisBCG);

  cout << "id = " << id << endl;

}
