import sys
import pysam


with pysam.VariantFile(sys.argv[1], threads=16) as f, pysam.VariantFile(sys.argv[2], 'w', threads=16, header=f.header) as f_out:
    for rec in f.fetch():
        len_ref = len(rec.alleles[0])
        len_alt = len(rec.alleles[1])
        if max(len_ref, len_alt) >= 50:
            if 'SVTYPE' in rec.info:
                f_out.write(rec)
            else:
                new_rec = rec.copy()
                new_rec.info['REFLEN'] = len_ref
                if len_ref == 1 and len_alt > 1:
                    new_rec.info['SVTYPE'] = 'INS'
                    new_rec.info['SVLEN'] = len_alt - len_ref
                elif len_ref > 1 and len_alt == 1:
                    new_rec.info['SVTYPE'] = 'DEL'
                    new_rec.info['SVLEN'] = len_alt - len_ref
                else:
                    new_rec.info['SVTYPE'] = 'COMPLEX'
                    new_rec.info['SVLEN'] = len_ref
                f_out.write(new_rec)
        else:
            f_out.write(rec)
        
