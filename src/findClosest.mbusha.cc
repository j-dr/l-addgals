#include "galaxy.h"
#include "particle.h"
#include "galaxy.h"
#include "halo.h"
#include "color.h"
#include "stl_util.h"
#include "myrand.h"
//#include "choose.h"
//#include "iostream.h"
#include <iostream>
#include<ctime>
#include<time.h>
#include <string.h>

extern double normal_random(float mean, float stddev);

class keyValue
{
public:
	int key;
	float value;
	keyValue(int key, float value)
	{
		this->key = key;
		this->value = value;
	}
	keyValue()
	{

	}	
};

class intkeyValue
{
public:
        int key;
        int value;
        intkeyValue(int key, float value)
        {
                this->key = key;
                this->value = value;
        }
        intkeyValue()
        {

        }
};

void *Mix(int *tab1,int *tab2,int count1,int count2)
//void Mix(int *tab1,int *tab2,int count1,int count2)
{
  int i,i1,i2;
  i = i1 = i2 = 0;
  int* temp = (int*)malloc(sizeof(int)*(count1+count2));

  while((i1<count1) && (i2<count2))
  {
    while((i1<count1) && ( (*(tab1+i1)) <= (*(tab2+i2))))
    {
      *(temp+i++) = *(tab1+i1);
      i1++;
    }
    if (i1<count1)
    {
      while((i2<count2) && ((*(tab2+i2)) <= (*(tab1+i1))))
      {
        *(temp+i++) = *(tab2+i2);
        i2++;
      }
    }
  }

  memcpy(temp+i,tab1+i1,(count1-i1)*sizeof(int));
  memcpy(tab1,temp,count1*sizeof(int));

  memcpy(temp+i,tab2+i2,(count2-i2)*sizeof(int));
  memcpy(tab2,temp+count1,count2*sizeof(int));

  free(temp);
}

void MergeSort(int* tab,int count)
{
  if (count==1) return;

  MergeSort(tab,count/2);
  MergeSort(tab+count/2,(count+1)/2);
  Mix(tab,tab+count/2,count/2,(count+1)/2);
}

void Mix(keyValue *tab1,keyValue *tab2,int count1,int count2)
{
  int i,i1,i2;
  i = i1 = i2 = 0;
  keyValue* temp = (keyValue *)malloc(sizeof(keyValue)*(count1+count2));
 
  while((i1<count1) && (i2<count2))
  {
    while((i1<count1) && ( (*(tab1+i1)).value <= (*(tab2+i2)).value ))
    {
      *(temp+i++) = *(tab1+i1);
      i1++;
    }
    if (i1<count1)
    {
      while((i2<count2) && ((*(tab2+i2)).value <= (*(tab1+i1)).value ))
      {
        *(temp+i++) = *(tab2+i2);
        i2++;
      }
    }
  }
 
  memcpy(temp+i,tab1+i1,(count1-i1)*sizeof(keyValue));
  memcpy(tab1,temp,count1*sizeof(keyValue));
 
  memcpy(temp+i,tab2+i2,(count2-i2)*sizeof(keyValue));
  memcpy(tab2,temp+count1,count2*sizeof(keyValue));
  
  free(temp);
}
 
void MergeSort(keyValue *tab,int count)
{
  if (count==1) return;
 
  MergeSort(tab,count/2);
  MergeSort(tab+count/2,(count+1)/2);
  Mix(tab,tab+count/2,count/2,(count+1)/2);
}

void Mix(intkeyValue *tab1,intkeyValue *tab2,int count1,int count2)
{
  int i,i1,i2;
  i = i1 = i2 = 0;
  intkeyValue* temp = (intkeyValue *)malloc(sizeof(intkeyValue)*(count1+count2));

  while((i1<count1) && (i2<count2))
  {
    while((i1<count1) && ( (*(tab1+i1)).value <= (*(tab2+i2)).value ))
    {
      *(temp+i++) = *(tab1+i1);
      i1++;
    }
    if (i1<count1)
    {
      while((i2<count2) && ((*(tab2+i2)).value <= (*(tab1+i1)).value ))
      {
        *(temp+i++) = *(tab2+i2);
        i2++;
      }
    }
  }

  memcpy(temp+i,tab1+i1,(count1-i1)*sizeof(intkeyValue));
  memcpy(tab1,temp,count1*sizeof(intkeyValue));

  memcpy(temp+i,tab2+i2,(count2-i2)*sizeof(intkeyValue));
  memcpy(tab2,temp+count1,count2*sizeof(intkeyValue));

  free(temp);
}

void MergeSort(intkeyValue *tab,int count)
{
  if (count==1) return;

  MergeSort(tab,count/2);
  MergeSort(tab+count/2,(count+1)/2);
  Mix(tab,tab+count/2,count/2,(count+1)/2);
}

void print(keyValue *tab, int count)
{
	for(int i=0;i<count;i++)
	{
		cout<<tab[i].key<<" "<<tab[i].value<<endl;
	}	
}


int binarySearch(keyValue* sortedArray, int first, int last, float key)
{	
	while (first <= last) 
	{
       		int mid = (first + last) / 2;  // compute mid point.
	       	if (key > sortedArray[mid].value) 
           		first = mid + 1;  // repeat search in top half.
       		else if (key < sortedArray[mid].value) 
           		last = mid - 1; // repeat search in bottom half.
       		else
           		return mid;     // found it. return position /////
   	}
	return first;
}

int binarySearch(keyValue* sortedArray, int last, float key)
{	
	return binarySearch(sortedArray, 0,  last,  key);
}


int binarySearch(float* sortedArray, int first, int last, float key)
{	
	while (first <= last) 
	{
       		int mid = (first + last) / 2;  // compute mid point.
	       	if (key > sortedArray[mid]) 
           		first = mid + 1;  // repeat search in top half.
       		else if (key < sortedArray[mid]) 
           		last = mid - 1; // repeat search in bottom half.
       		else
           		return mid;     // found it. return position /////
   	}
	return first;
}

void insert(int* temp1, int position, int key)
{
	int i = position-1;
	while(temp1[i] > key && i>=0)
	{
		temp1[i+1] = temp1[i];
		i--;
	}
	temp1[i+1] = key;
}

int findCloseGalaxies2(vector <GalSED> &v, float mag, float dens, float ThisZ, int ThisBCG)
{
	vector<GalSED>::iterator begin = v.begin();
	vector<GalSED>::iterator end = v.end();
	//#define INSERTIONSORT
	//#define ARRAYSIZE 10000
	#define ARRAYSIZE 20000

#ifdef RED_FRACTION
	static int STEP=1500; //what I have been using for BCC
	//static int STEP=6000;  ///testing to see how it impacts the blue cloud width
#else
	static int STEP=500;
#endif
	static int sorted = 0;
	static keyValue* densities;
	static keyValue* magnitudes;
	static int size;
	static int temp1[ARRAYSIZE];
	static int temp2[ARRAYSIZE];


	int answer;

	//cout<<"Searching for SED for galaxy with mag = "<<mag<<" dens = "<<dens<<endl;

#ifdef RED_FRACTION
	float emag = evolve_mag(mag, ThisZ);
	if (emag < -22) emag = -22;
	if (emag > -18) emag = -18;
	float mv = emag + 20.0;
	float zv = 1 / (ThisZ + 1) - 0.47;
	float red_fraction = fq0 + fqz1 * zv + fqz2 * zv*zv + fqz3 * zv*zv*zv + 
	                       fq1 * mv + fq2 * mv*mv + fq3 * mv*mv*mv +
	                       fq1z1 * mv*zv + fq1z2 * mv*zv*zv + fq2z1*mv*mv*zv;
	if (red_fraction > 1.0) red_fraction = 1.0;
	if (red_fraction < 0.0) red_fraction = 0.0;
#endif

	if( sorted != 1)
	{
		sorted = 1;
		size = end - begin;
		densities = new keyValue[size];
		magnitudes = new keyValue[size];
		vector<GalSED>::iterator inputIterator;
		int i;
		for ( i=0,inputIterator= begin ; inputIterator < end; inputIterator++, i++ )
    		{
			densities[i].key = i;
			magnitudes[i].key = i;
			densities[i].value = log10((*inputIterator).Dens());
			magnitudes[i].value = (*inputIterator).MR();
		}
		//cout<<"Starting magnitude sort"<<endl;
		MergeSort(magnitudes,size);
		//cout<<"Starting density sort"<<endl;
		MergeSort(densities,size);
		//cout<<"Done with both sorts"<<endl;
	}	
	int magIndex = binarySearch(magnitudes,size,mag);
	int densIndex = binarySearch(densities,size,log10(dens));

	//cout<<"Initial magIndex: "<<magIndex<<" value: "<<magnitudes[magIndex].value<<endl;
	//cout<<"Initial densIndex: "<<densIndex<<" value: "<<(pow(10.0, densities[densIndex].value))<<endl;

	int maxSteps = (int)(size/STEP);

	int count1=0;
	int count2=0;
	int stepIndex;

	insert(temp1, count1++, magnitudes[magIndex].key);
	insert(temp2, count2++, densities[densIndex].key);

	for(int m=0;m<maxSteps;m++)
	{
		for(int n=1;n<=STEP;n++)
		{
			if(count1>=ARRAYSIZE || count2>=ARRAYSIZE)
			{
				cout<<"count1= "<<count1<<"count2= "<<count2<<endl;
				answer = -1; 
				return answer;
			}
			stepIndex = m*STEP + n;
			#ifdef INSERTIONSORT
				if(magIndex+stepIndex < size)
					insert(temp1, count1++, magnitudes[magIndex+stepIndex].key);
				if(magIndex-stepIndex > 0)
					insert(temp1, count1++, magnitudes[magIndex-stepIndex].key);
				if(densIndex+stepIndex < size)
					insert(temp2, count2++, densities[densIndex+stepIndex].key);
				if(densIndex-stepIndex > 0)
					insert(temp2, count2++, densities[densIndex-stepIndex].key);
			#else
                                if(magIndex+stepIndex < size)
                                        temp1[count1++] = magnitudes[magIndex+stepIndex].key;
                                if(magIndex-stepIndex > 0)
                                        temp1[count1++] = magnitudes[magIndex-stepIndex].key;
                                if(densIndex+stepIndex < size)
                                        temp2[count2++] = densities[densIndex+stepIndex].key;
                                if(densIndex-stepIndex > 0)
                                        temp2[count2++] = densities[densIndex-stepIndex].key;
			#endif
		}

		#ifndef INSERTIONSORT
			MergeSort(temp1, count1);
			MergeSort(temp2, count2);
		#endif
		
                int i=0,j=0,k=0;
#ifdef RED_FRACTION
		//cout<<"Counting the number of matching red and blue galaxies..."<<endl;
		float nred = 0, ntot = 0;
                while(i < count1 && j < count2)
                {
                        if ( temp1[i] == temp2[j])
                        {
				ntot++;
				if(v[temp1[i]].Red()) nred++;
                                i++;
                                j++;
                        }
                        else if(temp1[i] > temp2[j])
                        {
                                j++;
                        }
                        else if(temp1[i] < temp2[j])
                        {
                                i++;
                        }
                }
		//cout<<" Total number of matching objects: "<<ntot<<", total number of red objects: "<<nred<<endl;
		float local_red_fraction = nred/ntot;
		//float target_local_red_fraction = local_red_fraction*red_fraction/REDFRACTION1;
		float target_local_red_fraction = local_red_fraction*red_fraction;

#ifdef RF_TEST	   
		string filename = "rftest.dat";
		ofstream rffile(filename.c_str(), std::ofstream::out | std::ofstream::app);

		if (rffile.fail()) {
		  cerr<<"Error:  cannont open Redfraction file file: "<<filename<<endl;
		  exit(1);
		}

		rffile << red_fraction << " " << local_red_fraction << " " << target_local_red_fraction << " " << ntot << endl;
#endif
	        int is_red = 1;
        	float ran = drand48();
        	if (ran > target_local_red_fraction)
                	is_red = 0;
        	//cout<<"   local_red_fraction = "<<local_red_fraction<<", target_local_red_fraction = "<<target_local_red_fraction<<", this is_red = "<<is_red<<endl;
		i=0,j=0,k=0;
#endif

		int answerCount = 0;
		int answerArray[ARRAYSIZE];
                while(i < count1 && j < count2)
                {
                        if ( temp1[i] == temp2[j])
                        {
#ifdef RED_FRACTION
                                if(v[temp1[i]].Red() == is_red)
                                {
                                        answer = temp1[i];
					answerArray[answerCount++] = answer;
                                }
#else
                                answer = temp1[i];
				answerArray[answerCount++] = answer;
#endif
                                i++;
                                j++;
                                k++;
                        }
                        else if(temp1[i] > temp2[j])
                        {
                                j++;
                        }
                        else if(temp1[i] < temp2[j])
                        {
                                i++;
                        }
                }
		if (answerCount > 0)
		{
			int ind = rand() % answerCount;
			return answerArray[ind];
		}
	}
	cout<<"Error!  Couldn't find a SED Match for galaxy with mag = "<<mag<<" dens = "<<dens<<endl;
	return -1;
}



//void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies, ChainEl chel, int LLBins, int * &NInBin, int * &BinStart)
#ifdef SHAM_TEST
void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos)
#else
void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies)
#endif
{


	//while(1)
	//{

	//int n; //dummy

	//This value is used for the number of galaxies that should be considered based on distance. Analogous to the size of the bin. Higher value results in lesser processing time.
	int NEARESTZ = 500000;
	float dZ = 0.01; //v2.11
	//float dZ = 1.005; //trying to improve v4.20
	//float dZ = 0.015; //MGS
#ifdef SNAPSHOT
	dZ = 2.0;  //We dont' care about redshift if we're dealing with a box
#endif
	
	//Thus value is for the maximum number of galaxies that will be considered. Generally, 100 is a good number. Higher value takes more processing time.
	//int MAXSEARCHD = galaxies.size()/100;  //for v2.11
	int MAXSEARCHD = galaxies.size()/1000;  //trying to improve v4.20
	float MINSEARCHD = 0.1;
	//int MAXSEARCHD = particles.size()/50; //for MGS
	//float MAXDELTA8 = 0.1;
	//int MAXSEARCHD = 5000;
	cout<<"MAXSEARCHD: "<<MAXSEARCHD<<endl;

	//How many buckets/bins do we need based on Z
	int BUCKETSZ = 10; 

	//cout<<"Press 0 to stop assigning galaxies."<<endl;
	//cin>>n;
	//if(n == 0)
	//	break;
	//cout<<"Number of Buckets?"<<endl;
	//cin>>BUCKETSZ;
	//cout<<"MAXSEARCHD?"<<endl;
	//cin>>MAXSEARCHD;
	//cout<<"NEARESTZ?"<<endl;
	//cin>>NEARESTZ;
	
	float* bucketVals = new float[BUCKETSZ];
	float* minBucketVals = new float[BUCKETSZ];
	float* maxBucketVals = new float[BUCKETSZ];

	keyValue* galaxyZ = new keyValue[galaxies.size()];
	keyValue* galaxyD = new keyValue[galaxies.size()];
	keyValue* galaxyMr = new keyValue[galaxies.size()]; //mbusha added to loop by magnitude
	int galaxyCount = galaxies.size();
	for(int i=0;i<galaxyCount;i++)
	{
		galaxyZ[i].key = i;
		galaxyZ[i].value = galaxies[i]->zGal();
		galaxyD[i].key = i;
		galaxyD[i].value = galaxies[i]->Dist8();
	}
	MergeSort(galaxyD,galaxyCount);
	//this lets us sort through in magnitude, but still use the sorted density list
	for(int i=0;i<galaxyCount;i++)
	{
		galaxyMr[i].key = i;
		galaxyMr[i].value = galaxies[galaxyD[i].key]->Mr();
	}
	MergeSort(galaxyMr,galaxyCount);
	
	keyValue* particleZ = new keyValue[particles.size()];
	keyValue* particleD = new keyValue[particles.size()];
	keyValue* particleZSorted = new keyValue[particles.size()];
	int particleCount = particles.size();
	for(int i=0;i<particleCount;i++)
	{
		particleZ[i].key = i;
		particleZ[i].value = particles[i]->Zred();
		particleZSorted[i].key = i;
		particleZSorted[i].value = particles[i]->Zred();
		particleD[i].key = i;
		particleD[i].value = particles[i]->Dist8();
	}
	MergeSort(particleZSorted,particleCount);
	MergeSort(particleD,particleCount);
	
	cout<<"ParticleCount= "<<particleCount<<" GalaxyCount= "<<galaxyCount<<endl;
	cout<<"First few sorted particle densities: "<<particleD[0].value<<" "<<particleD[1].value<<" "<<particleD[2].value<<endl;
	
	int bucketLength = (int)particleCount/BUCKETSZ;
	if(NEARESTZ < bucketLength)
	{
		NEARESTZ = bucketLength;
	}
	for(int i=0;i<BUCKETSZ;i++)
	{
		bucketVals[i] = particleZSorted[bucketLength*i].value;
		if(bucketLength*i - NEARESTZ >= 0 && i>0)
			minBucketVals[i] = particleZSorted[bucketLength*i - NEARESTZ].value;
		else
			minBucketVals[i] = particleZSorted[0].value;
		if(bucketLength*i + NEARESTZ < particleCount && i<BUCKETSZ-1)
			maxBucketVals[i] = particleZSorted[bucketLength*i + NEARESTZ].value;
		else
			maxBucketVals[i] = particleZSorted[particleCount-1].value;
	}
	
	bool* assigned = new bool[2*particleCount];
	//vector<bool> assigned;
	int NPreviouslyAssigned = 0;
	for(int i=0;i<particleCount;i++)
	{
		//assigned.push_back(false);
		assigned[i] = false;
		///*
		//if (particles[particleD[i].key]->IsGal()){
		if (particles[i]->IsGal()){
		  assigned[i] = true;
		  NPreviouslyAssigned++;
		}
		//*/
	}
	/*
	for(int i=0;i<galaxyCount;i++)
	  {
	    if(galaxies[i]->Central())
	      {
		assigned[galaxies[i]->P()->Pid()] = true;
		NPreviouslyAssigned++;
	      }
	  }
	*/
	cout<<"# of particles previously assigned to centrals: "<<NPreviouslyAssigned<<endl;

	int position=0;
	int gID;
	float gD;
	float gZ;
	//int bucketNumber;
	float minZVal;
	float maxZVal;
	int particleID;
	int pID;
	int pi;
	bool found;


	float percent = 0.0;
	int Multiple_assigned = 0;
	int NRandom = 0;
	int NRandomZ = 0;
	int NRandomZBright = 0;
	for(int gi=0;gi<galaxyCount;gi++)
	//for(int gi=0;gi<galaxyCount-0.1*galaxyCount;gi++)
	//for(int gi=0;gi<100;gi++)
	{
		//cout<<gi<<endl;
		if(((float)gi)/galaxyCount > percent)
		{
			cout<<100*percent<<"% done"<<endl;
			//percent += 0.01;
			percent += 0.1;
		}


		///mbusha question:  Why are we looping ghrough based on density?  Why high-dens first?
			//a:  makes for a faster Binary Search (line ~590).  Changed to loop over mag, 
			//    but preserver density sorting
		//gID = galaxyD[gi].key;
		//gD = galaxyD[gi].value;
                gID = galaxyD[galaxyMr[gi].key].key;
                gD = galaxyD[galaxyMr[gi].key].value;
		gZ = galaxyZ[gID].value;
		//if(gi < 100) cout<<galaxies[gID]->Mr()<<endl;

#ifdef SHAM_TEST
                //do we keep the initial galaxy or reassign it?
                if (read_hod == 1 && gD < rnn_cut && 
		    galaxies[gID]->Mhost() > mhost_cut){
                        //cout<<"Skipping galaxy..."<<endl;
                        //save the galaxy location and velocity to the particle
			Halo * thalo = galaxies[gID]->H();
                        Point xx(thalo->X()/(sim.LengthUnit()),
                                 thalo->Y()/(sim.LengthUnit()),
                                 thalo->Z()/(sim.LengthUnit()));
                        Point vv(thalo->Vx(),
                                 thalo->Vy(),
                                 thalo->Vz());
                        //Point xx(halos[gID]->X()/(sim.LengthUnit()),
                        //         halos[gID]->Y()/(sim.LengthUnit()),
                        //         halos[gID]->Z()/(sim.LengthUnit()));
                        //Point vv(halos[gID]->Vx(),
                        //         halos[gID]->Vy(),
                        //         halos[gID]->Vz());
			//cout<<"Created points..."<<endl;
			//xx.Print();
			//vv.Print();
                        Particle * particle = new Particle(xx,vv,gD);
			//cout<<"Created particle..."<<endl;
			//if (halos[gID]->Host() < 0) 
			if (thalo->Host() == thalo->Id()) 
			  galaxies[gID]->DefineCentral();

                        //add the particle and assign the galaxy to it
                        particles.push_back(particle);
			//cout<<"Added particle to vector..."<<endl;
                        particleID = particles.size() - 1;
			//cout<<"Updating assigned array -- expect segfault..."<<endl;
                        //assigned.push_back(true);
                        assigned[particleID] = true;
			//cout<<"Linking galaxy to particle..."<<endl;
                        galaxies[gID]->P(particles[particleID]);
			//cout<<"Defining particle as a galaxy..."<<endl;
                        particles[particleID]->MakeGal(gID);
			//cout<<"Fake particle created and assigned."<<endl;
			//if (gID < 10){
			//	xx.Print();
			//	galaxies[gID]->Print();
			//	particles[particleID]->PosPrint();
			//}
                        continue;
                }
#endif


		//mbusha change:  Skip an object if it has been assigned as a central
		if (galaxies[gID]->Central())
		  continue;

/*
		//This piece of code is not being used when the other piece of code is being used
		bucketNumber = binarySearch(bucketVals,0,BUCKETSZ,gZ);
		minZVal = minBucketVals[bucketNumber];
		maxZVal = maxBucketVals[bucketNumber];

		//Alternative code, takes more time than the above, but is expected to yield better results
		bucketNumber = binarySearch(particleZSorted,0,particleCount,gZ);
		if(bucketNumber-NEARESTZ > 0)
			minZVal = particleZSorted[bucketNumber-NEARESTZ].value;
		else
			minZVal = particleZSorted[0].value;
		if(bucketNumber+NEARESTZ < particleCount)
			maxZVal = particleZSorted[bucketNumber+NEARESTZ].value;
		else
			maxZVal = particleZSorted[particleCount-1].value;
*/

		//Michael's version -- just define an explicit dZ
		minZVal = gZ - dZ;
		maxZVal = gZ + dZ;

		//original version -- requires that we loop in order of increasing density
		//position = binarySearch(particleD, position, particleCount,gD);
		//mbusha's new version -- lets up look via magnitude instead
		position = binarySearch(particleD, 0, particleCount,gD);

		//cout<<"MinZVal="<<minZVal<<" MaxZVal="<<maxZVal<<endl;
		particleID = particleD[position].key;
		found = false;
		pi = -1;
		float diff1 = 0.0;
		float diff2 = 0.0;
		while(pi<MAXSEARCHD-1 && diff1 < MINSEARCHD && diff2 < MINSEARCHD)
		//for(pi=0;pi<MAXSEARCHD;pi++)
		  //for(pi=0;pi<particles.size();pi++)
		{
			pi++;
			if( position-pi >= 0)
			{
				pID = particleD[position-pi].key;
				if (pID <= 0 || pID >= particleCount) continue;
                                if (read_hod)
                                  if (particles[pID]->Mhost() > mhost_cut && particles[pID]->Dist8() < rnn_cut) 
                                    continue;
				diff1 = fabs(gD-particleD[position-pi].value);
				//if(fabs(particleD[position].value - particleD[position-pi].value) > MAXDELTA8) break;
				if(!assigned[pID] && particleZ[pID].value>=minZVal && particleZ[pID].value<=maxZVal)
				{
					particleID = pID;
					found = true;
					break;
				}
			}

			if( position+pi < particleCount-1)
			{
				pID = particleD[position+pi].key;
				if (pID <= 0 || pID >= particleCount) continue;
                                if (read_hod)
                                  if (particles[pID]->Mhost() > mhost_cut && particles[pID]->Dist8() < rnn_cut) 
                                    continue;
				diff2 = fabs(gD-particleD[position-pi].value);
				//if(fabs(particleD[position].value - particleD[position-pi].value) > MAXDELTA8) break;
				if(!assigned[pID] && particleZ[pID].value>=minZVal && particleZ[pID].value<=maxZVal)
				{
					found = true;
					particleID = pID;
					break;
				}
			}
		}
		
		//cout<<"GalaxyZ="<<galaxies[gID]->zGal()<<" GalaxyD="<<galaxies[gID]->Dist8()<<" ParticleZ="<<particles[particleID]->Zred()<<" ParticleD="<<particles[particleID]->Dist8()<<endl;
		//if(found)
			//cout<<"Found at "<<pi<<endl;
		//else
			//cout<<"Not found.......................... :( "<<endl;
		if(!found){
		  if (galaxies[gID]->Mr() < -18.0){
		    cout<<"  ERROR!  Particle not found for galaxy "<<gi<<" with d8 = "<<gD<<", z = "<<gZ<<", Mr = "<<galaxies[gID]->Mr()<<endl;
		    cout<<"     zmin = "<<minZVal<<", zmax = "<<maxZVal<<endl;
		    //cout<<"     d8_min = "<<particleD[position-MAXSEARCHD].value<<", d8_max = "<<particleD[position+MAXSEARCHD].value<<", d8_desired = "<<particleD[position].value<<endl;
		  }
		  for(pi=0;pi<particles.size();pi++)
		    {
		      if( position-pi >= 0)
			{
			  pID = particleD[position-pi].key;
			  if (pID <= 0 || pID >= particleCount) continue;
			  if (read_hod) 	
			    if (particles[pID]->Mhost() > mhost_cut && particles[pID]->Dist8() < rnn_cut) 
			      continue;
			  if(!assigned[pID] && particles[pID]->Save())
			    {
			      particleID = pID;
			      found = true;
			      break;
			    }
			}
		      
		      if( position+pi < particleCount)
			{
			  pID = particleD[position+pi].key;
			  if (pID <= 0 || pID >= particleCount) continue;
                          if (read_hod)
                            if (particles[pID]->Mhost() > mhost_cut && particles[pID]->Dist8() < rnn_cut) 
                              continue;
			  if(!assigned[pID] && particles[pID]->Save())
			    {
			      found = true;
			      particleID = pID;
			      break;
			    }
			}
		    }
		  if(!found){
		    cout<<" Couldn't find any decent match. Randomly Assigning."<<endl;
		    NRandom++;
		  }
		  if (galaxies[gID]->Mr() < -18.0){
		    cout<<"  going to assign particle with d8 = "<<particles[particleID]->Dist8()<<", z = "<<particles[particleID]->Zred()<<" (gi = "<<gi<<")"<<endl;
		    cout<<"  ``randomly'' selected particles current assignment state: "<<assigned[particleID]<<endl;
		    NRandomZBright++;
		  }
		  NRandomZ++;
		  if(assigned[particleID]) Multiple_assigned++;
		}
		assigned[particleID] = true;
		galaxies[gID]->P(particles[particleID]);
		particles[particleID]->MakeGal(gID);
	}

	cout<<"Finished assignment.  "<<endl;
	cout<<endl;
	cout<<Multiple_assigned<<" galaxies had duplicate particles (out of "<<galaxies.size()<<" galaxies)."<<endl;
	cout<<NRandomZ<<" were randomly assigned in z. "<<endl;
	cout<<NRandomZBright<<"  bright galaxies were randomly assigned in z. "<<endl;
	cout<<NRandom<<" were just plain shitty matches."<<endl;
	cout<<"First Few galaxy positions: "<<endl;
	cout<<"   "<<galaxies[0]->P()->X()<<" "<<galaxies[0]->P()->Y()<<" "<<galaxies[0]->P()->Z()<<endl;
	cout<<"   "<<galaxies[1]->P()->X()<<" "<<galaxies[1]->P()->Y()<<" "<<galaxies[1]->P()->Z()<<endl;
	cout<<"   "<<galaxies[2]->P()->X()<<" "<<galaxies[2]->P()->Y()<<" "<<galaxies[2]->P()->Z()<<endl;
	cout<<"   "<<galaxies[3]->P()->X()<<" "<<galaxies[3]->P()->Y()<<" "<<galaxies[3]->P()->Z()<<endl;
	cout<<"   "<<galaxies[4]->P()->X()<<" "<<galaxies[4]->P()->Y()<<" "<<galaxies[4]->P()->Z()<<endl;
	cout<<endl;
}

void Read_L_BCG(float &M0, float &Mc, float &a, float &b, float &k)
{
  string filename = lbcgfile;
/*
  string pdf_base = lbcgfile;
  string filename =  "denspdf/"+pdf_base+"099.txt";
  if (ZREDMIN > 0.027)
    filename = "denspdf/"+pdf_base+"098.txt";
  if (ZREDMIN > 0.054)
    filename = "denspdf/"+pdf_base+"097.txt";
  if (ZREDMIN > 0.11)
    filename = "denspdf/"+pdf_base+"095.txt";
  if (ZREDMIN > 0.2)
    filename = "denspdf/"+pdf_base+"092.txt";
  if (ZREDMIN > 0.27)
    filename = "denspdf/"+pdf_base+"090.txt";
  if (ZREDMIN > 0.33)
    filename = "denspdf/"+pdf_base+"088.txt";
  if (ZREDMIN > 0.44)
    filename = "denspdf/"+pdf_base+"085.txt";
  if (ZREDMIN > 0.56)
    filename = "denspdf/"+pdf_base+"082.txt";
  if (ZREDMIN > 0.64)
    filename = "denspdf/"+pdf_base+"080.txt";
  if (ZREDMIN > 0.73)
    filename = "denspdf/"+pdf_base+"078.txt";
  if (ZREDMIN > 0.87)
    filename = "denspdf/"+pdf_base+"075.txt";
  if (ZREDMIN > 1.14)
    filename = "denspdf/"+pdf_base+"070.txt";
  if (ZREDMIN > 1.25)
    filename = "denspdf/"+pdf_base+"068.txt";
  if (ZREDMIN > 1.43)
    filename = "denspdf/"+pdf_base+"065.txt";
  if (ZREDMIN > 1.77)
    filename = "denspdf/"+pdf_base+"060.txt";
*/

  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"Error:  cannont open LBCG file: "<<filename<<endl;
    exit(1);
  }
  else
    cout<<"Reading LBCG file "<<filename<<endl;

  char line[200];
  //file>>line;
  //  file.getline(line, 200);
  cout<<"Line: "<<endl;
  file>>M0>>Mc>>a>>b>>k;
  cout<<M0<<" "<<Mc<<" "<<a<<" "<<b<<" "<<k<<endl;
  file.close();
}

void LinkHalosParticles(vector <Particle *> &P, vector <Halo *> &H)
{
  cout<<"Linking Halos with their Most Bound Particles..."<<endl;
  intkeyValue* hids = new intkeyValue[H.size()];
  intkeyValue* pids = new intkeyValue[P.size()];

  cout<<"Setting up P Keys..."<<endl;
  //setup our lookup tables for particle id #'s in poth particle and halo arrays
  for(int i=0;i<P.size();i++){
    pids[i].key = i;
    pids[i].value = P[i]->Pid();
  }
  cout<<"Setting up H Keys..."<<endl;
  for(int i=0;i<H.size();i++){
    hids[i].key = i;
    hids[i].value = H[i]->Particle();
  }

  cout<<"Sorting the tables..."<<endl;
  //sort our lookup tables
  MergeSort(pids,P.size());
  MergeSort(hids,H.size());

  //do a single loop through lookup table to find particles at halo centers
  cout<<"Looping through particles..."<<endl;
  int this_ih = 0;
  for(int ip=0;ip<P.size();ip++)
    {
      //we found a match -- assign particle and move to the next halo
      //cout<<ip<<endl;
      if (pids[ip].value == hids[this_ih].value)
	{
	  //cout<<"Found a match for halo "<<hids[this_ih].key<<" to particle id "<<pids[ip].value<<"...";
	  H[hids[this_ih].key]->Particle(pids[ip].key);
	  P[pids[ip].key]->Hid(hids[this_ih].key);
	  //cout<<"assigned the particle."<<endl;
	  //cout<<"  Halo Position info:  ";
	  //H[hids[this_ih].key]->Print();
	  //cout<<"  Particle Position info:  ";
	  //P[pids[ip].key]->PosPrint();
	  this_ih++;
	  if (this_ih == H.size())
	    break;
	}
      //this is an error -- particle isn't in our list!!!
      if (pids[ip].value > hids[this_ih].value)
	{
	  //if the halo is outside our volume we make a new particle
	  //cout<<"Couldn't find a particle for halo "<<this_ih<<endl;
	  int ih = hids[this_ih].key;
	  if(H[ih]->InVol())
	    {
	      //cout<<"Error!  Didn't readin particle "<<hids[this_ih].value<<" that should be the center for halo "<<hids[this_ih].key<<endl;
	      //cout<<"Just making up a new particle to go here."<<endl;
	    }
	  float xfac = 1./(sim.LengthUnit());
	  Point xx(H[ih]->X()*xfac,H[ih]->Y()*xfac,H[ih]->Z()*xfac);
	  Point vv(H[ih]->Vx(),H[ih]->Vy(),H[ih]->Vz());
	  float dist8 = H[ih]->Dist8();
	  Particle * particle = new Particle(xx,vv,dist8);
	  particle->MakeGal(99999999);
	  particle->Hid(ih);
	  particle->MVir(H[ih]->M());
	  particle->RVir(H[ih]->R200());
	  particle->RHalo(0);
	  P.push_back(particle);
	  H[ih]->Particle(P.size()-1);
	  this_ih++;
	  if (this_ih == H.size())
	    break;
	  ip--; //make sure we continue searching at this particle in our loop
	  //if(H[ih]->InVol())
	  //{
	  //cout<<"Particle added.  this_ih = "<<this_ih<<" of "<<H.size()<<endl;
	  //}
	}
    }
  cout<<"Checking for extra halos to add..."<<endl;
  //Did we not cover some of the halos?  If so, add new particles.  
  if(this_ih < H.size()){
    cout<<"Making up particles for halos "<<this_ih<<" to "<<H.size()-1<<endl;
    int start = this_ih;
    for(this_ih=start;this_ih<H.size();this_ih++){
      int ih = hids[this_ih].key;
      float xfac = 1./(sim.LengthUnit());
      Point xx(H[ih]->X()*xfac,H[ih]->Y()*xfac,H[ih]->Z()*xfac);
      Point vv(H[ih]->Vx(),H[ih]->Vy(),H[ih]->Vz());
      float dist8 = H[ih]->Dist8();
      Particle * particle = new Particle(xx,vv,dist8);
      particle->Hid(ih);
      particle->MVir(H[ih]->M());
      particle->RVir(H[ih]->R200());
      particle->RHalo(0);
      particle->MakeGal(99999999);
      P.push_back(particle);
      H[ih]->Particle(P.size()-1);
    }
  }
  cout<<"Finished Linking Halo Particles."<<endl;

}

double CalculateMeanBCGLum(double lnLc0, double ALc, double Mpiv, double BLc, double MVir, double z)
{

  return lnLc0 + ALc * log(MVir / Mpiv) + BLc * log(1 + z);

}

void AssignBCGs(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos)
{
	//int n; //dummy

	//This value is used for the number of galaxies that should be considered based on magnitude. Analogous to the size of the bin. Higher value results in lesser processing time.
	int NEARESTMR = 5000;
#ifdef SNAPSHOT
	NEARESTMR = 25000;
#endif
	float dMr = 0.15;
	//float dMr = 0.05;

	//Thus value is for the maximum number of galaxies that will be considered. Generally, 100 is a good number. Higher value takes more processing time.
	//int MAXSEARCHD = 5000;
	int MAXSEARCHD = galaxies.size()/100;
	//int MAXSEARCHD = galaxies.size()/30;
	//int MAXSEARCHD = galaxies.size(); //find closest density with magnitude difference within dMr
#ifdef SHAM_TEST
	MAXSEARCHD = 0;
#endif

	//How many buckets/bins do we need based on Z
	int BUCKETSMR = 1000; 
	//int BUCKETSMR = galaxies.size()/NEARESTMR; 

	cout<<"Assigning BCGs with NEARESTMR = "<<NEARESTMR<<", MAXSEARCHD = "<<MAXSEARCHD<<", BUCKETSMR = "<<BUCKETSMR<<endl;

	string outcheck="./BCG_check.dat";
	ofstream bcg_dens_file(outcheck.c_str());

	//cout<<"Press 0 to stop assigning galaxies."<<endl;
	//cin>>n;
	//if(n == 0)
	//	break;
	//cout<<"Number of Buckets?"<<endl;
	//cin>>BUCKETSZ;
	//cout<<"MAXSEARCHD?"<<endl;
	//cin>>MAXSEARCHD;
	//cout<<"MEARESTMR?"<<endl;
	//cin>>MEARESTMR;

	LinkHalosParticles(particles, halos);
	
	//set magnitudes based on halo mass
	cout<<"Assigning BCG magnitudes..."<<endl;

#ifdef DESCLF
        double sigma_L = 0.364;
        double lnLc0   = 24.554;
        double ALc     = 0.355;
        double BLc     = 0.936;
        double Mpiv    = 2.35*1e14;
	double Msunr   = 4.67;
#else
	float M0, Mc, a, b, k;
	Read_L_BCG(M0, Mc, a, b, k);
        Mc = pow(10, Mc);
#endif 

	for(int i=0;i<halos.size();i++)
	  {
#ifdef SHAM_TEST
	    if (halos[i]->Host() != halos[i]->Id() || halos[i]->M() < BCG_Mass_lim) {
		halos[i]->Mr(99);
		continue;
            }
#else
	    //	    if (halos[i]->Host() >= 0 || halos[i]->M() < BCG_Mass_lim) {
	    if (halos[i]->M() < BCG_Mass_lim) {	    
		halos[i]->Mr(99);
		//continue;
            }
#endif
	    //	    if (halos[i]->Host() < 0 && halos[i]->M() >= BCG_Mass_lim) {	
	    if (halos[i]->M() >= BCG_Mass_lim) {    
#ifdef DESCLF
	      //assign using power law fit to DES CLF
	      double lnL0 = CalculateMeanBCGLum(lnLc0, ALc, Mpiv, BLc, halos[i]->M(), halos[i]->ZredReal());
	      //add scatter
	      double lnL  = normal_random(lnL0, sigma_L);
	      double mr = -2.5 * log10(exp(lnL)) + Msunr;
#else
	      double m200 = halos[i]->M();
	      if (m200 > 1e15) m200 = 1e15;
	      double mr0 = M0 - 2.5*(a*log10(m200/Mc) - b*log10(1.+pow(m200/Mc,k/b)));
	      //float scatter = 0.17;
	      double mr = normal_random(mr0, 2.5*SCATTER);
#endif
	      halos[i]->Mr(mr);
	      
	    }
#ifdef DEBUG
	    if(i<10) cout<<"Halo "<<i<<": M200 = "<<m200<<", Mr = "<<mr<<endl;
#endif

	    /*
	    //assign through my SHAM fit
	    double fit0 = fit_00 + fit_01*halos[hi]->Zred();
	    double fit1 = fit_10 + fit_11*halos[hi]->Zred();
	    mr = fit0 + fit1*log10(m200);
	    lum = pow(10.0, -0.4*(mr - Mstar));
	    lum_before = lum;
	    lum = pow(10.0, normal_random(log10(lum), 0.15));
	    mr = -2.5*log10(lum) + Mstar;
	    //halos[i]->Mr(mr);
	    */
	  }

	//float* bucketVals = new float[BUCKETSMR];
	//float* minBucketVals = new float[BUCKETSMR];
	//float* maxBucketVals = new float[BUCKETSMR];

	cout<<"Generating Halo Keys..."<<endl;
	keyValue* haloMr = new keyValue[halos.size()];
	keyValue* haloD = new keyValue[halos.size()];
	int haloCount = halos.size();
	for(int i=0;i<haloCount;i++)
	{
		haloMr[i].key = i;
		haloMr[i].value = halos[i]->Mr();
		haloD[i].key = i;
		haloD[i].value = halos[i]->Dist8();
	}
	MergeSort(haloD,haloCount);
	
	cout<<"Generating Galaxy Keys..."<<endl;
	keyValue* galaxyMr = new keyValue[galaxies.size()];
	keyValue* galaxyD = new keyValue[galaxies.size()];
	//keyValue* galaxyMrSorted = new keyValue[galaxies.size()];
	int galaxyCount = galaxies.size();
	for(int i=0;i<galaxyCount;i++)
	{
		galaxyMr[i].key = i;
		galaxyMr[i].value = galaxies[i]->Mr();
		//galaxyMrSorted[i].key = i;
		//galaxyMrSorted[i].value = galaxies[i]->Mr();
		galaxyD[i].key = i;
		galaxyD[i].value = galaxies[i]->Dist8();
	}
	//cout<<"Sorting Magnitudes..."<<endl;
	//MergeSort(galaxyMrSorted,galaxyCount);
	cout<<"Sorting Densities... ("<<galaxyD[0].value<<", "<<galaxyD[1].value<<", "<<galaxyD[2].value<<"...)"<<endl;
/*
	float maxD = -10.;
	float minD = 10.;
	int maxD_ind = -1;
	int minD_ind = -1;
	ofstream fout("particle_densities.txt");
	for(int i=0;i<galaxyCount;i++)
	{
		fout<<galaxyD[i].value<<"\n";
		if (galaxyD[i].value > maxD)
		{
			maxD = galaxyD[i].value;
			maxD_ind = galaxyD[i].key;
		}
		if (galaxyD[i].value < minD)
		{
			minD = galaxyD[i].value;
			minD_ind = galaxyD[i].key;
		}
	}
	fout.close();
	cout<<"Min/Max galaxy density and indicies: "<<minD<<" ("<<minD_ind<<"), "<<maxD<<" ("<<maxD_ind<<")"<<endl;
*/
	MergeSort(galaxyD,galaxyCount);
	
	cout<<"galaxyCount= "<<galaxyCount<<" HaloCount= "<<haloCount<<endl;
        cout<<" some sorted galaxy densities: "<<galaxyD[0].value<<" "<<galaxyD[1].value<<" "<<galaxyD[2].value<<"... "<<galaxyD[galaxyCount-1].value<<endl;
	
	/*
	int bucketLength = (int)galaxyCount/BUCKETSMR;
	if(NEARESTMR < bucketLength)
	{
		NEARESTMR = bucketLength;
	}
	for(int i=0;i<BUCKETSMR;i++)
	{
		bucketVals[i] = galaxyMrSorted[bucketLength*i].value;
		if(bucketLength*i - NEARESTMR >= 0 && i>0)
			minBucketVals[i] = galaxyMrSorted[bucketLength*i - NEARESTMR].value;
		else
			minBucketVals[i] = galaxyMrSorted[0].value;
		if(bucketLength*i + NEARESTMR < galaxyCount && i<BUCKETSMR-1)
			maxBucketVals[i] = galaxyMrSorted[bucketLength*i + NEARESTMR].value;
		else
			maxBucketVals[i] = galaxyMrSorted[galaxyCount-1].value;
	}
	*/

	bool* assigned = new bool[galaxyCount];
	for(int i=0;i<galaxyCount;i++)
	{
		assigned[i] = false;
	}

	int position=0;
	int hID;
	float hD;
	float hMr;
	//int bucketNumber;
	float minMrVal;
	float maxMrVal;
	int galaxyID;
	int gID;
	int gi;
	bool found;
	int ExtraGals = 0;
	int hInVol = 0;

	float percent = 0.0;

	//int ibad = 1882048;
	int ibad = 2082048000;
	for(int hi=0;hi<haloCount;hi++)
	{
	        if (hi >= ibad) 
		{
			cout<<"Halo of "<<hi<<" of "<<haloCount<<endl;
			cout<<"x/y/z = "<<halos[hi]->X()<<" "<<halos[hi]->Y()<<" "<<halos[hi]->Z()<<endl;
			cout<<"ra/dec = "<<halos[hi]->Ra()<<" "<<halos[hi]->Dec()<<endl;
		}
		if(((float)hi)/haloCount > percent)
		{
			cout<<100*percent<<"% done"<<endl;
			//percent += 0.01;
			percent += 0.1;
		}

		hID = haloD[hi].key;
		hD = haloD[hi].value;
		hMr = haloMr[hID].value;
		//if (hMr >= 0) continue;

		if (hi >= ibad) cout<<"Have halo properties.  Checking location."<<endl;
		if(!(halos[hID]->InVol()))
		{
		  if (hi >= ibad) cout<<"Halo is NOT in the volume."<<endl;
		  continue;
		}
		hInVol++;
		if (hi >= ibad) cout<<"Halo is in the volume."<<endl;

		/*
		//This piece of code is not being used when the other piece of code is being used
		bucketNumber = binarySearch(bucketVals,0,BUCKETSMR,hMr);
		minMrVal = minBucketVals[bucketNumber];
		maxMrVal = maxBucketVals[bucketNumber];

		//Alternative code, takes more time than the above, but is expected to yield better results
		bucketNumber = binarySearch(galaxyMrSorted,0,galaxyCount,hMr);
		if(bucketNumber-NEARESTMR > 0)
			minMrVal = galaxyMrSorted[bucketNumber-NEARESTMR].value;
		else
			minMrVal = galaxyMrSorted[0].value;
		if(bucketNumber+NEARESTMR < galaxyCount)
			maxMrVal = galaxyMrSorted[bucketNumber+NEARESTMR].value;
		else
			maxMrVal = galaxyMrSorted[galaxyCount-1].value;
		*/

		position = binarySearch(galaxyD, position, galaxyCount,hD);
		//cout<<"  Halo rnn: "<<hD<<", closest galaxy rnn: "<<galaxyD[position].value<<endl;

		//mbusha's alternate:  set min and max by specified boundaries
		minMrVal = hMr - dMr;
		maxMrVal = hMr + dMr;

		//cout<<"MinMrVal="<<minMrVal<<" MaxMrVal="<<maxMrVal<<", position = "<<position<<endl;
		galaxyID = galaxyD[position].key;
		found = false;
		//cout<<"Searching galaxiees..."<<endl;
		if (hi >= ibad) cout<<"Searching for galaxies..."<<endl;
		for(gi=0;gi<MAXSEARCHD;gi++)
		{
			if( position-gi >= 0)
			{
				gID = galaxyD[position-gi].key;
				if(!assigned[gID] && galaxyMr[gID].value>=minMrVal && galaxyMr[gID].value<=maxMrVal)
				{
					galaxyID = gID;
					found = true;
					break;
				}
			}

			if( position+gi < galaxyCount)
			{
				gID = galaxyD[position+gi].key;
				if(!assigned[gID] && galaxyMr[gID].value>=minMrVal && galaxyMr[gID].value<=maxMrVal)
				{
					found = true;
					galaxyID = gID;
					break;
				}
			}
		}
		//cout<<"Finished search.  found = "<<found<<endl;
		//cout<<"HaloM="<<halos[hID]->Mr()<<" HaloD="<<halos[gID]->Dist8()<<" GalaxyM="<<galaxies[gID]->Mr()<<" GalaxyD="<<galaxies[gID]->Dist8()<<", found = "<<found<<endl;
		//if(found)
		//cout<<"Found at "<<gID<<endl;
		//else
		//cout<<"Not found.......................... :( "<<endl;
		if(found)
		  {
		    if (hi >= ibad) cout<<"Found galaxy with ID "<<galaxyID<<endl;
		    assigned[galaxyID] = true;
		  }
		else
		  {
		    //create a new galaxy for the halo center
		    if (hi >= ibad) cout<<"Creating new galaxy."<<endl;
		    Galaxy * galaxy = new Galaxy(hMr,galaxies.size(),hD);
		    galaxy->Dist8(halos[hID]->Dist8());
		    galaxies.push_back(galaxy);
		    galaxyID = galaxies.size()-1;
		    ExtraGals++;
		  }
		//halos[hID]->P(particles[particleID]);
		//galaxies[galaxyID]->MakeGal(gID);

		//cout<<"Setting halo/galaxy info..."<<endl;
		if (hi >= ibad) cout<<"Setting halo/galaxy info..."<<endl;
		int pid = halos[hID]->Particle();
		galaxies[galaxyID]->P(particles[pid]);
		particles[pid]->MakeGal(galaxyID);
		galaxies[galaxyID]->DefineCentral();
		galaxies[galaxyID]->Mr(halos[hID]->Mr());
		
		//write out to our log file to see how well we did
		if (hi >= ibad) cout<<"Saving to log file..."<<endl;
		bcg_dens_file<<halos[hID]->Dist8()<<" "<<halos[hID]->Mr()<<" "<<halos[hID]->Zred()<<" "<<galaxies[galaxyID]->Dist8()<<" "<<galaxies[galaxyID]->Mr()<<" "<<galaxies[galaxyID]->zGal()<<" "<<found<<" "<<galaxyID<<endl;


	}

	cout<<"Did BCG assignment.  There were "<<hInVol<<" halos in the volume.  Had to create "<<ExtraGals<<" new galaxies ("<<((float) ExtraGals)*100.0/hInVol<<"%)"<<endl;
	//}
}
