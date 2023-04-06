from subprocess import call

if __name__ == '__main__':
    PS1 = [1500. + 100 * x for x in range(51)]
    PS2 = [3500.]
    den_scale = [1.4]
    te_scale = [2.]
    angle = [-35.]
    for ii in PS1:
        for jj in PS2:
            for kk in den_scale:
                for ll in te_scale:
                    for mm in angle:
                        fread = open('genray_original.in', 'r')
                        fwrite = open('genray.in', 'w')
                        for line in fread:
                            array = line.split()
                            if len(array) > 0:
                                if array[0] == 'curc=':
                                    array[1] = str(jj*40.)
                                    array[2] = str(ii*40.)
                                    array[3] = str(ii*40.)
                                    array[4] = str(jj*40.)
                                    array[5] = str(jj*40.)
                                    array[6] = str(jj*40.)
                                    fwrite.write(" ".join(array))
                                    fwrite.write("\n")
                                elif array[0] == 'temp_scale=':
                                    array[1] = str(ll)
                                    fwrite.write(" ".join(array))
                                    fwrite.write("\n")
                                elif array[0] == 'den_scale=':
                                    array[1] = str(kk)
                                    fwrite.write(" ".join(array))
                                    fwrite.write("\n")
                                elif array[0] == 'betast=':
                                    array[1] = str(mm)
                                    fwrite.write(" ".join(array))
                                    fwrite.write("\n")
                                else:
                                    fwrite.write(line)
                            else:
                                fwrite.write(line)
                        fread.close()
                        fwrite.close()
                        log = open('log', 'w')
                        call(['/home/icl/genray-c/genray-c_160826.2/xgenray'], shell = False, stdout = log)
                        ncfile = "_".join(['genray',str(ii),str(jj),str(kk), str(ll),str(mm),'real_col.nc'])
                        call(['mv', 'genray.nc', ncfile], shell = False)
    

    
