#include "loadVCF.h"

int getVCFHeaderAll(const char *fname)
{
    int i;
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) error("Could not read %s\n", fname);
    
    tbx_t *tbx = tbx_index_load(fname);
    if ( !tbx ) error("Could not load .tbi/.csi index of %s\n", fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(fp, KS_SEP_LINE, &str) >= 0 )
        {
        if ( !str.l || str.s[0]!=tbx->conf.meta_char ) break;
        puts(str.s);
    }
    free(str.s);
    tbx_destroy(tbx);
    if ( hts_close(fp) ) error("hts_close returned non-zero status: %s\n", fname);    
    return 0;
}


int getVCFHeader(const char *fname)
{
    int i;
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) error("Could not read %s\n", fname);
    
    tbx_t *tbx = tbx_index_load(fname);
    if ( !tbx ) error("Could not load .tbi/.csi index of %s\n", fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(fp, KS_SEP_LINE, &str) >= 0 )
        {
        if ( !str.l || str.s[0]!=tbx->conf.meta_char ) break;
        if(str.s[1]!='#') puts(str.s);
    }
    free(str.s);
    tbx_destroy(tbx);
    if ( hts_close(fp) ) error("hts_close returned non-zero status: %s\n", fname);
    return 0;
}

/*
int main(int argc, char** argv){
    getVCFHeader(argv[1]);
    return 0;
}
*/
