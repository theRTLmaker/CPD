#include "physics.h"
#include <stdio.h>




particle_t* addParticle(particle_t *head, particle_t *p){
    p -> nextParticle = head;
    head = p;
    return head;
}




particle_t * removeParticle(particle_t *head, particle_t *p){
    particle_t *cur;
    particle_t *aux;

    if(p == head){
        head = p -> nextParticle;
        p -> nextParticle = NULL;
        return head;
    }

    aux = head;
    cur = head -> nextParticle;

    while(cur != NULL){
        if(cur == p){
            aux -> nextParticle = p -> nextParticle;
            p -> nextParticle = NULL;
            break;
        }
        aux = aux-> nextParticle;
        cur = cur -> nextParticle;
    }

    return head;
}

