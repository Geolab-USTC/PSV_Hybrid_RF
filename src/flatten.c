#define EARTH_RADIUS     6371.0 
double flatten(int layer, float *tth)
{
    int i;
    double thickness = 0.0;

    for(i=0; i<layer-1; i++)
        thickness += tth[i];
    thickness += 0.5*tth[layer-1]; 

    return (EARTH_RADIUS/(EARTH_RADIUS-thickness));
}


