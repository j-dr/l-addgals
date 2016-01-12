#ifndef findclosest_h
#define findclosest_h

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

void MergeSort(keyValue *tab,int count);

#endif
