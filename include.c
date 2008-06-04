#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*****************************************************************************/
int treat_line(char *line) {

   char *ptr;
   
   if (strncmp(line + 6, "INCLUDE", 6) != 0) {
      printf("%s\n", line);
      return(0);
   }

   ptr = rindex(line, '\'');
   if (ptr) *ptr = 0;
   ptr = index(line, '\'');
   if (ptr) treat_file(ptr + 1);
   
   
   return(0);
}

/*****************************************************************************/
int treat_file(char *filename) {
   char line[1024];
   FILE *fp;
   
   /* Open file */
   
   fp = fopen(filename, "r");
   if (!fp) {
      fprintf(stderr, "Unable to read file %s - %m\n", filename);
      return(-1);
   }
   
   /* Read file */
   
   while(1) {
      if (fgets(line, sizeof(line) - 1, fp) == NULL) break;
      line[strlen(line) - 1] = 0;
      treat_line(line);
   }
   
   /* Close file */
   
   fclose(fp);
   
   return(0);
}

/*****************************************************************************/
int main(int argc, char **argv) {

   int i;
   for (i = 1; i < argc; i++) treat_file(argv[i]);
   return(0);
}
