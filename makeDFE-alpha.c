/*
 Version 2020.01.04
 
 Copyright (C) 2020 Tuomas Hamala
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 For any other inquiries send an email to tuomas.hamala@gmail.com
 
 Program for producing SFS and divergence-counts reguired by DFE-alpha
 
 compiling: gcc makeDFE-alpha.c -o makeDFE-alpha -lm
 
 usage:
 -coord [file] coordinates file produced by 'show-coords' program from MUMmer (use settings -H -T)
 -div [file] substitution file produced by 'show-snps' program from MUMmer (use settings -C -I -H -T)
 -sites [file] tab-delimited file with chromosome and postition (0-fold or 4-fold)
 -vcf [file] full vcf-file containing variant and invariant sites
 -region [file] tab-delimited file with chromosome, start, and end for regions to use (optional)
 -gc [int] different DAF-classes 1 [WS] 2 [SW] 3 [SS] 4 [WW] 5 [SS+WW]
 
 example:
 ./makeDFE-alpha -coord thaliana-lyrata.filt.coord -div thaliana-lyrata.filt.snps -sites 4fold.sites -vcf thaliana.full.vcf -gc 1 > out.4fold.WS.txt
 
 The program was written for A.thaliana data. It assumes that all files are sorted by chromosome and position, chromosomes are identified with numbers (e.g. 1 and not chr1), and VCF-file contains no heterozygote sites.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#define merror "\nERROR: System out of memory\n"

typedef struct{
    int chr, start, stop;
}Region_s;

typedef struct{
    int chr, pos;
    char ref, alt;
}Site_s;

void openFiles(int argc, char *argv[]);
Region_s *readCoord(FILE *coord_file, int *n);
Region_s *readTarget(FILE *target_file, int *n);
int **readSites(FILE *site_file, Region_s *coords, Region_s *target, int coord_n, int target_n, int *n);
Site_s *readDiv(FILE *div_file, int **sites, int site_n, int *n);
void readVcf(FILE *vcf_file, int **sites, Site_s *div, int gc, int site_n, int div_n);
void lineTerminator(char *line);

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
    timer = time(NULL);
    
    openFiles(argc, argv);
    
    second = time(NULL) - timer;
    
    minute = second / 60;
    
    hour = second / 3600;
    
    if(isatty(1))
        fprintf(stderr,"\n");
    
    if(hour > 0)
        fprintf(stderr,"Run finished in %i h, %i min & %i sec\n\n", hour, minute-hour*60, second-minute*60);
    else if(minute > 0)
        fprintf(stderr,"Run finished in %i min & %i sec\n\n", minute, second-minute*60);
    else if(second > 5)
        fprintf(stderr,"Run finished in %i sec\n\n", second);
    else
        fprintf(stderr,"\n");
    
    return 0;
}

void openFiles(int argc, char *argv[]){
    
    int i, gc=0, coord_n=0, div_n=0, site_n=0, target_n=0, **sites;
    char **list;
    FILE *vcf_file=NULL, *coord_file=NULL, *div_file=NULL, *site_file=NULL, *target_file=NULL;
    Region_s *coords, *target;
    Site_s *div;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-coord") == 0){
            if((coord_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-coord %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-div") == 0){
            if((div_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-div %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-sites") == 0){
            if((site_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sites %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-vcf") == 0){
            if((vcf_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-vcf %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-region") == 0){
            if((target_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-region %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-gc") == 0){
            gc = atoi(argv[++i]);
            if(gc < 0 || gc > 5){
                fprintf(stderr,"\nERROR: allowed values for -gc are 1 [WS], 2 [SW] 3 [SS] 4 [WW] 5 [SS+WW]\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-gc %s\n", argv[i]);
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    
    fprintf(stderr,"\n");
    
    if(coord_file == NULL || div_file == NULL || site_file == NULL || vcf_file == NULL){
        fprintf(stderr,"ERROR: The following parameters are required: -coord [file] -div [file] -sites [file] -vcf [file]\n\n");
        exit(EXIT_FAILURE);
    }
    
    coords = readCoord(coord_file, &coord_n);
    if(target_file != NULL)
        target = readTarget(target_file, &target_n);
    sites = readSites(site_file, coords, target, coord_n, target_n, &site_n);
    div = readDiv(div_file, sites, site_n, &div_n);
    readVcf(vcf_file, sites, div, gc, site_n, div_n);
}

Region_s *readCoord(FILE *coord_file, int *n){
    
    int i, line_i=0;
    char c, *line=NULL, *temp=NULL;
    Region_s *list;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(coord_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(coord_file);
    
    if((list = malloc(line_i*sizeof(Region_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, coord_file)) != -1){
        lineTerminator(line);
        list[*n].start = atoi(strtok(line,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        for(i=0;i<6;i++)
            temp = strtok(NULL,"\t");
        if(isdigit(temp[0])){
            list[*n].chr = atoi(temp);
            *n = *n + 1;
        }
    }
    
    free(line);
    fclose(coord_file);
    
    return list;
}

Region_s *readTarget(FILE *target_file, int *n){
    
    int i, line_i=0;
    char c, *line=NULL;
    Region_s *list;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(target_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(target_file);
    
    if((list = malloc(line_i*sizeof(Region_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, target_file)) != -1){
        lineTerminator(line);
        list[*n].chr = atoi(strtok(line,"\t"));
        list[*n].start = atoi(strtok(NULL,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        *n = *n + 1;
    }
    
    free(line);
    fclose(target_file);
    
    return list;
}

int **readSites(FILE *site_file, Region_s *coords, Region_s *target, int coord_n, int target_n, int *n){
    
    int i, j, line_i=0, brk=0, coord_i=0, target_i=0, **list;
    char c, *line=NULL;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(site_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(site_file);
    
    if((list = malloc(line_i*sizeof(int*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<line_i;i++){
        if((list[i] = malloc(2*sizeof(int))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while((read = getline(&line, &len, site_file)) != -1){
        lineTerminator(line);
        if(isdigit(line[0]) == 0)
            continue;
        list[*n][0] = atoi(strtok(line,"\t"));
        list[*n][1] = atoi(strtok(NULL,"\t"));
        while(coord_i < coord_n){
            if(list[*n][0] == coords[coord_i].chr){
                if(list[*n][1] >= coords[coord_i].start && list[*n][1] <= coords[coord_i].stop){
                    if(target_n > 0){
                        brk = 0;
                        while(target_i < target_n){
                            if(list[*n][0] == target[target_i].chr){
                                if(list[*n][1] >= target[target_i].start && list[*n][1] <= target[target_i].stop){
                                    *n = *n + 1;
                                    brk = 1;
                                    break;
                                }
                                else if(list[*n][1] < target[target_i].start)
                                    break;
                            }
                            else if(list[*n][0] < target[target_i].chr)
                                break;
                            target_i++;
                        }
                        if(brk == 1)
                            break;
                    }
                    else{
                        *n = *n + 1;
                        break;
                    }
                }
                else if(list[*n][1] < coords[coord_i].start)
                    break;
            }
            else if(list[*n][0] < coords[coord_i].chr)
                break;
            coord_i++;
        }
    }
    
    free(coords);
    free(line);
    fclose(site_file);
    
    return list;
}

Site_s *readDiv(FILE *div_file, int **sites, int site_n, int *n){
    
    int i, line_i=0, site_i=0;
    char c, *line=NULL, *temp=NULL;
    Site_s *list;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(div_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(div_file);
    
    if((list = malloc(line_i*sizeof(Site_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, div_file)) != -1){
        lineTerminator(line);
        list[*n].pos = atoi(strtok(line,"\t"));
        temp = strtok(NULL,"\t");
        list[*n].ref = temp[0];
        temp = strtok(NULL,"\t");
        list[*n].alt = temp[0];
        for(i=0;i<6;i++)
            temp = strtok(NULL,"\t");
        if(isdigit(temp[0])){
            list[*n].chr = atoi(temp);
            while(site_i < site_n){
                if(list[*n].chr == sites[site_i][0]){
                    if(list[*n].pos == sites[site_i][1]){
                        *n = *n + 1;
                        break;
                    }
                    else if(list[*n].pos < sites[site_i][1])
                        break;
                }
                else if(list[*n].chr < sites[site_i][0])
                    break;
                site_i++;
            }
        }
    }
    
    free(line);
    fclose(div_file);
    
    return list;
}

void readVcf(FILE *vcf_file, int **sites, Site_s *div, int gc, int site_n, int div_n){
    
    int i, j, chr=0, pos=0, i0=0, i1=0, count=0, ind_i=0, site_i=0, div_i=0, ok=0, rd=0, *geno;
    double di=0, si=0, *sfs;
    char ref, alt, *line=NULL, *temp=NULL;
    size_t len=0;
    ssize_t read;
    
    while((read = getline(&line, &len, vcf_file)) != -1){
        lineTerminator(line);
        
        temp = strtok(line,"\t");
        
        if(strcmp(temp, "#CHROM") == 0){
            i = 1;
            while(temp != NULL){
                if(i > 9)
                    ind_i++;
                temp = strtok(NULL,"\t");
                i++;
            }
            
            if((geno=malloc(ind_i*sizeof(int))) == NULL){
                fprintf(stderr,merror);
                exit(EXIT_FAILURE);
            }
            
            if((sfs=malloc((ind_i+1)*sizeof(double))) == NULL){
                fprintf(stderr,merror);
                exit(EXIT_FAILURE);
            }
            
            for(i=0;i<=ind_i;i++)
                sfs[i] = 0;
        }
        else if(isdigit(temp[0])){
            chr = atoi(temp);
            i = 0;
            j = 1;
            i1 = 0;
            i0 = 0;
            ok = 0;
            rd = 0;
            while(temp != NULL){
                if(j == 2){
                    pos =  atoi(temp);
                    while(site_i < site_n){
                        if(chr == sites[site_i][0]){
                            if(pos == sites[site_i][1]){
                                while(div_i < div_n){
                                    if(chr == div[div_i].chr){
                                        if(pos == div[div_i].pos){
                                            rd = 1;
                                            break;
                                        }
                                        else if(pos < div[div_i].pos)
                                            break;
                                    }
                                    else if(chr < div[div_i].chr)
                                        break;
                                    div_i++;
                                }
                                ok = 1;
                                break;
                            }
                            else if(pos < sites[site_i][1])
                                break;
                        }
                        else if(chr < sites[site_i][0])
                            break;
                        site_i++;
                    }
                    if(ok == 0)
                        break;
                }
                else if(j == 4){
                    ref = temp[0];
                    if(rd == 1 && ref != div[div_i].ref){
                        fprintf(stderr,"Warning: ref alleles differ at chr %i pos %i\n", chr, pos);
                        ok = 0;
                        break;
                    }
                }
                else if(j == 5){
                    alt = temp[0];
                    if(gc != 0){
                        ok = 0;
                        if(rd == 1 && alt != div[div_i].alt)
                            break;
                        if(ref == '.' && alt == '.')
                            break;
                        if(gc == 3){
                            if((ref == 'G' || ref == 'C' || ref == '.') && (alt == 'G' || alt == 'C' || alt == '.'))
                                ok = 1;
                        }
                        else if(gc == 4){
                            if((ref == 'A' || ref == 'T' || ref == '.') && (alt == 'A' || alt == 'T' || alt == '.'))
                                ok = 1;
                        }
                        else if(gc == 5){
                            if((ref == 'G' || ref == 'C' || ref == '.') && (alt == 'G' || alt == 'C' || alt == '.'))
                                ok = 1;
                            else if((ref == 'A' || ref == 'T' || ref == '.') && (alt == 'A' || alt == 'T' || alt == '.'))
                                ok = 1;
                        }
                        else if(rd == 0){
                            if(gc == 1 && (ref == 'A' || ref == 'T' || ref == '.') && (alt == 'G' || alt == 'C' || alt == '.'))
                                ok = 1;
                            else if(gc == 2 && (ref == 'G' || ref == 'C' || ref == '.') && (alt == 'A' || alt == 'T' || alt == '.'))
                                ok = 1;
                        }
                        else if(rd == 1){
                            if(gc == 1 && (ref == 'G' || ref == 'C' || ref == '.') && (alt == 'A' || alt == 'T' || alt == '.'))
                                ok = 1;
                            else if(gc == 2 && (ref == 'A' || ref == 'T' || ref == '.') && (alt == 'G' || alt == 'C' || alt == '.'))
                                ok = 1;
                        }
                        if(ok == 0)
                            break;
                    }
                }
                else if(j > 9){
                    if(temp[0] == '1' && temp[2] == '1'){
                        geno[i] = 1;
                        i1++;
                    }
                    else if(temp[0] == '0' && temp[2] == '0'){
                        geno[i] = 0;
                        i0++;
                    }
                    else
                        geno[i] = 9;
                    i++;
                }
                temp = strtok(NULL,"\t");
                j++;
            }
            if(ok == 0)
                continue;
            for(i=0;i<ind_i;i++){
                if(rd == 0 && geno[i] == 1)
                    count++;
                else if(rd == 0 && geno[i] == 9 && i1 > i0)
                    count++;
                else if(rd == 1 && geno[i] == 0)
                    count++;
                else if(rd == 1 && geno[i] == 9 && i0 > i1)
                    count++;
            }
            sfs[count]++;
            count = 0;
            if(rd == 1)
                di++;
            si++;
        }
    }
    
    for(i=0;i<=ind_i;i++)
        printf("%.0f ", sfs[i]);
    printf("\n");
    printf("%.0f %.0f\n", si, di);
    
    free(line);
    free(sfs);
    free(geno);
    fclose(vcf_file);
}

void lineTerminator(char *line){
    
    int i;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
