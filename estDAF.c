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
 
 Program for estimating derived allele frequencies
 
 compiling: gcc estDAF.c -o estDAF -lm
 
 usage:
 -coord [file] coordinates file produced by 'show-coords' program from MUMmer (use settings -H -T)
 -div [file] substitution file produced by 'show-snps' program from MUMmer (use settings -C -I -H -T)
 -vcf [file] vcf-file with variant sites
 -genes [file] tab-delimited file with name, chromosome, start, and end for each gene
 -gc [int] different DAF-classes 1 [WS] 2 [SW] 3 [SS] 4 [WW] 5 [SS+WW]
 
 example:
 ./estDAF -coord thaliana-lyrata.filt.coord -div thaliana-lyrata.filt.snps -vcf thaliana.poly.vcf -genes thaliana.genes.txt -gc 1 > out.WS.txt
 
 The program assumes that all files are sorted by chromosome and position, and that chromosomes are identified with numbers (e.g. 1 and not chr1)
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
    char id[50];
}Region_s;

typedef struct{
    int chr, pos;
    char ref, alt;
}Site_s;

void openFiles(int argc, char *argv[]);
Region_s *readCoord(FILE *coord_file, int *n);
Region_s *readGenes(FILE *gene_file, int *n);
Site_s *readDiv(FILE *div_file, int *n);
void readVcf(FILE *vcf_file, Region_s *genes, Region_s *coords, Site_s *div, int gc, int gene_n, int coord_n, int div_n);
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
    
    int i, gc=0, coord_n=0, div_n=0, site_n=0, vcf_n, gene_n=0;
    FILE *vcf_file=NULL, *coord_file=NULL, *div_file=NULL, *gene_file=NULL;
    Region_s *coords, *genes;
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
        
        else if(strcmp(argv[i], "-vcf") == 0){
            if((vcf_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-vcf %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-genes") == 0){
            if((gene_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-genes %s\n", argv[i]);
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
    
    if(gene_file == NULL || coord_file == NULL || div_file == NULL || vcf_file == NULL){
        fprintf(stderr,"ERROR: The following parameters are required: -coord [file] -div [file] -vcf [file] -genes [file]\n\n");
        exit(EXIT_FAILURE);
    }
    
    genes = readGenes(gene_file, &gene_n);
    coords = readCoord(coord_file, &coord_n);
    div = readDiv(div_file, &div_n);
    readVcf(vcf_file, genes, coords, div, gc, gene_n, coord_n, div_n);
}

Region_s *readGenes(FILE *gene_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0, coord_i=0;
    char c, *line=NULL;
    Region_s *list;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(gene_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(gene_file);
    
    if((list = malloc(line_i*sizeof(Region_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, gene_file)) != -1){
        lineTerminator(line);
        strncpy(list[*n].id, strtok(line,"\t"),49);
        list[*n].chr = atoi(strtok(NULL,"\t"));
        list[*n].start = atoi(strtok(NULL,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        *n = *n + 1;
    }
    
    free(line);
    fclose(gene_file);
    
    return list;
}

Region_s *readCoord(FILE *coord_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    char c, *line=NULL, *temp=NULL;
    Region_s *list;
    size_t len=0;
    ssize_t read;
    
    while((c=fgetc(coord_file)) != EOF){
        if(c == '\n')
            line_i++;
    }
    
    rewind(coord_file);
    
    if((list = malloc((line_i)*sizeof(Region_s))) == NULL){
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

Site_s *readDiv(FILE *div_file, int *n){
    
    int i, j, char_i=0, line_i=0, maxchar=0;
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
            *n = *n + 1;
        }
    }

    free(line);
    fclose(div_file);
    
    return list;
}

void readVcf(FILE *vcf_file, Region_s *genes, Region_s *coords, Site_s *div, int gc, int gene_n, int coord_n, int div_n){
    
    int i, chr=0, pos=0, gene_i=0, coord_i=0, div_i=0, ok=0, rd=0, da_i=0, a_i=0, s_i=0;
    double daf=0;
    char ref, alt, *line=NULL, *temp=NULL, *gene=NULL;
    size_t len=0;
    ssize_t read;
    
    while((read = getline(&line, &len, vcf_file)) != -1){
        lineTerminator(line);
        temp = strtok(line,"\t");
        if(isdigit(temp[0])){
            chr = atoi(temp);
            i = 1;
            while(temp != NULL){
                if(i == 2){
                    pos =  atoi(temp);
                    ok = 0;
                    rd = 0;
                    while(gene_i < gene_n){
                        if(chr == genes[gene_i].chr){
                            if(pos <= genes[gene_i].stop && pos >= genes[gene_i].start){
                                while(coord_i < coord_n){
                                    if(chr == coords[coord_i].chr){
                                        if(pos <= coords[coord_i].stop && pos >= coords[coord_i].start){
                                            if(gene == NULL){
                                                gene = genes[gene_i].id;
                                                printf("gene\tDAF\tnSites\n");
                                            }
                                            else if(strcmp(gene, genes[gene_i].id) != 0){
                                                daf = (double)da_i/(double)a_i;
                                                printf("%s\t%f\t%i\n", gene, daf, s_i);
                                                da_i = 0;
                                                a_i = 0;
                                                s_i = 0;
                                                gene = genes[gene_i].id;
                                            }
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
                                        else if(pos < coords[coord_i].start)
                                            break;
                                    }
                                    else if(chr < coords[coord_i].chr)
                                        break;
                                    coord_i++;
                                }
                                break;
                            }
                            else if(pos < genes[gene_i].start)
                                break;
                        }
                        else if(chr < genes[gene_i].chr)
                            break;
                        gene_i++;
                    }
                    if(ok == 0)
                        break;
                }
                else if(i == 4){
                    ref = temp[0];
                    if(rd == 1 && ref != div[div_i].ref){
                        fprintf(stderr,"Warning: ref alleles differ at chr %i pos %i\n", chr, pos);
                        break;
                    }
                }
                else if(i == 5){
                    alt = temp[0];
                    if(rd == 1 && alt != div[div_i].alt)
                        break;
                    if(gc != 0){
                        ok = 0;
                        if(gc == 3){
                            if((ref == 'G' && alt == 'C') || (ref == 'C' && alt == 'G'))
                                ok = 1;
                        }
                        else if(gc == 4){
                            if((ref == 'A' && alt == 'T') || (ref == 'T' && alt == 'A'))
                                ok = 1;
                        }
                        else if(gc == 5){
                            if((ref == 'G' && alt == 'C') || (ref == 'C' && alt == 'G'))
                                ok = 1;
                            else if((ref == 'A' && alt == 'T') || (ref == 'T' && alt == 'A'))
                                ok = 1;
                        }
                        else if(rd == 0){
                            if(gc == 1 && (ref == 'A' || ref == 'T') && (alt == 'G' || alt == 'C'))
                                ok = 1;
                            else if(gc == 2 && (ref == 'G' || ref == 'C') && (alt == 'A' || alt == 'T'))
                                ok = 1;
                        }
                        else{
                            if(gc == 1 && (ref == 'G' || ref == 'C') && (alt == 'A' || alt == 'T'))
                                ok = 1;
                            else if(gc == 2 && (ref == 'A' || ref == 'T') && (alt == 'G' || alt == 'C'))
                                ok = 1;
                        }
                        if(ok == 0)
                            break;
                    }
                }
                else if(i > 9){
                    if(rd == 1 && (temp[0] == '0' && temp[2] == '0'))
                        da_i++;
                    else if(rd == 0 && (temp[0] == '1' && temp[2] == '1'))
                        da_i++;
                    if(temp[0] != '.' && temp[2] != '.')
                        a_i++;
                }
                temp = strtok(NULL,"\t");
                i++;
                if(temp == NULL)
                    s_i++;
            }
        }
    }
    
    free(line);
    free(coords);
    free(genes);
    free(div);
    fclose(vcf_file);
}

void lineTerminator(char *line){
    
    int i;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
