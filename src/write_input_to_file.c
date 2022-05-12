//
// Created by Jitin Nair on 14/01/21.
//

#include "spinodal_omp.h"

void write_input_to_file(options_t options){

    //output label
    char* label = malloc( sizeof(char) * (strlen(options.label)+1) );
    strcpy( label, options.label);

    // create output directory if not created
    char path[60]="../output/";
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }

    strcat(path, label);
    strcat(path,"/");
    struct stat st2 = {0};
    if (stat(path, &st2) == -1) {
        mkdir(path, 0700);
    }

    // create folder within output directory named $label
    char input[30]="Input";
    char filename[strlen(path) + strlen(input) + 5];
    strcpy(filename, path);
    strcat(filename, input);
    strcat(filename, ".txt");

    FILE *fp;
    fp = fopen(filename,"w");

    fprintf(fp, "%-20s %s\n", "initial_conc", options.initial_conc);

    fclose(fp);

}