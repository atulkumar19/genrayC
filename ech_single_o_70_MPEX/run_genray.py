from subprocess import call

if __name__ == '__main__':
    B_scale = [1]
    den_scale = [1,.9,1.1]
    ax_launch = [-.153,-.123,-.093,-.063,-.033,-.003,.033,.063]
    rad_launch = [.03]
    phi_launch = [0.]
    angle = [-30,-20,-40]
    angle1 = [180]
    for ii in B_scale:
        for jj in den_scale:
            for kk in ax_launch:
                for ll in rad_launch:
                    for mm in phi_launch:
                        for nn in angle:
                            for oo in angle1:
                                fread = open('genray_original.in', 'r')
                                fwrite = open('genray.in', 'w')
                                for line in fread:
                                    array = line.split()
                                    if len(array) > 0:
                                        if array[0] == 'curc=':
                                            array[1] = str(ii*247080.3)
                                            array[2] = str(ii*348432.9)
                                            array[3] = str(ii*591373)
                                            array[4] = str(ii*429352.7)
                                            array[5] = str(ii*315939)
                                            array[6] = str(ii*344292.1)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")
                                        elif array[0] == 'zst=':
                                            array[1] = str(kk)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")
                                        elif array[0] == 'rst=':
                                            array[1] = str(ll)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")
                                        elif array[0] == 'phist=':
                                            array[1] = str(mm)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")    
                                        elif array[0] == 'den_scale=':
                                            array[1] = str(jj)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")
                                        elif array[0] == 'betast=':
                                            array[1] = str(nn)
                                            fwrite.write(" ".join(array))
                                            fwrite.write("\n")
                                        elif array[0] == 'alfast=':
                                            array[1] = str(oo)
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
                                ncfile = "_".join(['genray',str(ii),str(jj),str(kk), str(ll),str(mm),str(nn),str(oo),'real_col.nc'])
                                call(['mv', 'genray.nc', ncfile], shell = False)
    

    
