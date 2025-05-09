function bsave(fname, matrix, arg2)
% saves a binary file
% bsave(fname, matrix, arg2)
% uses command fwrite(fp, matrix, arg2);

fp = fopen(fname, 'w');

fwrite(fp, matrix, arg2);

fclose(fp);