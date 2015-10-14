#include "particle.h"
#include "galaxy.h"
#include "color.h"
#include "stl_util.h"
#include "myrand.h"
#include "choose.h"
#include "iostream.h"
#include<ctime>
#include<time.h>

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

void *Mix(int *tab1,int *tab2,int count1,int count2)
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
	binarySearch(sortedArray, 0,  last,  key);
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

int findCloseGalaxies2(vector <GalSED> &v, float mag, float dens)
{
	vector<GalSED>::iterator begin = v.begin();
	vector<GalSED>::iterator end = v.end();
	//#define INSERTIONSORT
	#define ARRAYSIZE 10000

	static int STEP=500;
	static int sorted = 0;
	static keyValue* densities;
	static keyValue* magnitudes;
	static int size;
	static int temp1[ARRAYSIZE];
	static int temp2[ARRAYSIZE];


	int answer;

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
		while(i < count1 && j < count2)	
		{
			if ( temp1[i] == temp2[j])
			{
				//cout<<"Equals! answer="<<temp1[i]<<" m= "<<m<<endl;
				answer = temp1[i];
				return answer;
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
	}
	return -1;
}



//void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies, ChainEl chel, int LLBins, int * &NInBin, int * &BinStart)
void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies)
{


	//while(1)
	//{

	int n; //dummy

	//This value is used for the number of galaxies that should be considered based on distance. Analogous to the size of the bin. Higher value results in lesser processing time.
	int NEARESTZ = 500000;
	
	//Thus value is for the maximum number of galaxies that will be considered. Generally, 100 is a good number. Higher value takes more processing time.
	int MAXSEARCHD = 1000;

	//How many buckets/bins do we need based on Z
	int BUCKETSZ = 10; 
#ifdef SNAPSHOT
	BUCKETSZ = 1;
#endif

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
	int galaxyCount = galaxies.size();
	for(int i=0;i<galaxyCount;i++)
	{
		galaxyZ[i].key = i;
		galaxyZ[i].value = galaxies[i]->zGal();
		galaxyD[i].key = i;
		galaxyD[i].value = galaxies[i]->Dist8();
	}
	MergeSort(galaxyD,galaxyCount);
	
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
	
	bool* assigned = new bool[particleCount];
	for(int i=0;i<particleCount;i++)
	{
		assigned[i] = false;
	}

	int position=0;
	int gID;
	float gD;
	float gZ;
	int bucketNumber;
	float minZVal;
	float maxZVal;
	int particleID;
	int pID;
	int pi;
	bool found;


	float percent = 0.0;

	for(int gi=0;gi<galaxyCount;gi++)
	{
		if(((float)gi)/galaxyCount > percent)
		{
			cout<<100*percent<<"% done"<<endl;
			//percent += 0.01;
			percent += 0.1;
		}

		gID = galaxyD[gi].key;
		gD = galaxyD[gi].value;
		gZ = galaxyZ[gID].value;

		//mbusha change:  Skip an object if it has been assigned as a central
		if (galaxies[gID]->Central())
		  continue;

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

		position = binarySearch(particleD, position, particleCount,gD);

		//cout<<"MinZVal="<<minZVal<<" MaxZVal="<<maxZVal<<endl;
		particleID = particleD[position].key;
		found = false;
		for(pi=0;pi<MAXSEARCHD;pi++)
		{
			if( position-pi >= 0)
			{
				pID = particleD[position-pi].key;
				if(!assigned[pID] && particleZ[pID].value>=minZVal && particleZ[pID].value<=maxZVal)
				{
					particleID = pID;
					found = true;
					break;
				}
			}

			if( position+pi < particleCount)
			{
				pID = particleD[position+pi].key;
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
		assigned[particleID] = true;
		galaxies[gID]->P(particles[particleID]);
		particles[particleID]->MakeGal(gID);
	}

	//}
}
