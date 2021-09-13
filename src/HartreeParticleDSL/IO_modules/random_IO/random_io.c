#include <stdlib.h>
#include "part.h"

struct part* random_io(int num_parts, struct config_type* config){
    struct part* parts = malloc(sizeof(struct part) * num_parts);
    for( int p=0; p < num_parts; p++){
        parts[p].core_part.position[0] = (double)(rand()) / (double)(RAND_MAX);
        parts[p].core_part.position[1] = (double)(rand()) / (double)(RAND_MAX);
        parts[p].core_part.position[2] = (double)(rand()) / (double)(RAND_MAX);
        parts[p].core_part.velocity[0] = 0.0;
        parts[p].core_part.velocity[1] = 0.0;
        parts[p].core_part.velocity[2] = 0.0;
    }
    config->space.nparts = num_parts;
    return parts;
}
