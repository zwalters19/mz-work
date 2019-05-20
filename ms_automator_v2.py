import os
import json



def get_spectra_matches():
    fls = os.listdir('.\\MAT3')
    output = ""
    for fl in fls:
        if fl.endswith("MAT"):
            output += fl.split('.')[0].split('_')[0] + '\t'
            mat_txt = ""
            with open(".\\MAT3\\"+fl,"r") as inpf:
                msms_flag = False
                skip_flag = False
                ms2_exp = ""
                for line in inpf:
                    if skip_flag == True:
                        words = line.split(" ")
                        ms2_exp += words[0] + "," + words[1][0:-1] + ";"
                    elif msms_flag == True:
                        skip_flag = True
                    #elif line[0:5] == "NAME:":
                        #i_d = line[5:-1].strip() + "\t"
                    elif line[0:11] == "MSTYPE: MS2":
                        msms_flag = True
                    elif line[0:13] == "PRECURSORMZ: ":
                        mz = line[13:-1] + "\t"
                    elif line[0:15] == "RETENTIONTIME: ":
                        rt = line[15:-1] + "\t"
                mat_txt += mz+rt+ms2_exp
            try:
                with open(".\\MAT3\\"+fl.split('.')[0]+"\\Spectral DB search.sfd","r") as s_file:
                    ms2_txt = ""
                    text = s_file.read()
                    if len(text) != 0:
                        should_split_next = False
                        entries = text.split('\n\n')
                        for line in entries[0].split('\n'):
                            if should_split_next == True:
                                if line != "\n":
                                    ms2_txt += line.split('\t')[0]+','+line.split('\t')[1]+';'
                                else:
                                    break
                            if line[:6]=="NAME: ":
                                ms2_txt += line[6:]+'\t'
                            if line[:12]=="TotalScore: ":
                                ms2_txt += line[12:]+"\t"
                            if line[:10]=="INCHIKEY: ":
                                ms2_txt += line[10:]+"\t"
                            if line[:5]=="Num P":
                                should_split_next = True
            except:
                ms2_txt = "\t\t\t;"
            if ms2_txt == "":
                ms2_txt = '\t\t\t;'
            output +=mat_txt[:-1]+"\t"+ms2_txt[:-1]+'\n'
        with open('spectra_matches.txt','w') as out_file:
            out_file.write(output)
    return 1
    
def msms_structures():
    count = 0
    output = list()
    output.append("ID\tM/Z\tRT\tMS2-Experimental\tAdduct\tMS2-Match\tName\tScore\tInchi\tFormula\tDMZ\tIso_score\tIso_match\tpossible_compounds\n")
    fl = os.listdir('.\\MAT2')
    
    f_c_matches = {}
    try:   
        file =  open('compound_dict.txt','r')
        f_c_matches = json.load(file)
        file.close()
    except:
        pass
    
    for file in fl:
        in_file = file.split('.')[0]
        i = in_file.split('_')[0]
        adduct = in_file.split('_')[1]
        if file.endswith('fgt'):
            best_line = ""
            best_flag = False
            msms_flag = False
            name = "\t"
            ms_score = "\t"
            score = "\t"
            ms2_exp = ""
            ms2_match = "\t"
            inchi = "\t"
            mz = "\t"
            rt = "\t"
            formula = '\t'
            dmz = "\t"
            isotope_match = "\t"
            with open(".\\MAT2\\"+in_file + ".MAT", "r") as inpf:
                msms_flag = False
                skip_flag = False
                for line in inpf:
                    if skip_flag == True:
                        words = line.split(" ")
                        ms2_exp += words[0] + "," + words[1][0:-1] + ";"
                    elif msms_flag == True:
                        skip_flag = True
                    elif line[0:11] == "MSTYPE: MS2":
                        msms_flag = True
                    elif line[0:13] == "PRECURSORMZ: ":
                        mz = line[13:-1] + "\t"
                    elif line[0:15] == "RETENTIONTIME: ":
                        rt = line[15:-1] + "\t"
            the_dir =  ".\\MAT2\\" +  in_file
            files = os.listdir(the_dir)
            if len(files) != 0:
                for j in files:
                    compound = ""
                    formula = j[:-4]
                    if formula in f_c_matches.keys():
                        compound = str(f_c_matches[formula])
                    with open('.\\MAT2\\'+in_file+".fgt",'r') as f_file:
                        file = f_file.read()
                        entries = file.split('\n\n')
                        matched = False
                        for entry in entries:
                            for line in entry.split('\n'):
                                if matched == True:
                                    if line [:11] == "TOTALSCORE:":
                                        score = line[12:-1] + "\t"
                                    if line[:15] == 'MASSDIFFERENCE:':
                                        dmz = line[16:-1]+"\t"
                                    if line[:24] == "ISOTOPICINTENSITY[M+1]: ":
                                        isotope_match = "1;" + line[24:]
                                    if line[:24] == "ISOTOPICINTENSITY[M+2]: ":
                                        isotope_match = isotope_match + ";" + line[24:]
                                        break
                                if line[:6] == "NAME: ":
                                    if line[6:] != formula:
                                        break
                                    else:
                                        matched = True
                    s_file = the_dir + "\\" + j
                    with open(s_file, 'r') as inpf:
                        s_file = inpf.read()
                        if len(s_file) > 0:
                            entries = s_file.split("\n\n")
                            for entry in entries:
                                ms2_match = ""
                                lines = entry.split("\n")
                                for line in lines:
                                    if msms_flag == True:
                                        words = line.split("\t")
                                        if len(words) > 2:
                                            ms2_match += words[2] + "," + str(100) + ";"
                                        else:
                                            pass
                                    elif line[0:6] == "NAME: ":
                                        name = line[6:] + '\t'
                                    elif line[0:10] == "INCHIKEY: ":
                                        inchi = line[10:] + '\t'
                                    elif line[0:12] == "TotalScore: ":
                                        ms_score = line[12:] + "\t"
                                    elif line[0:3] == "Num":
                                        msms_flag = True
                                msms_flag = False
                                if best_flag == True:
                                    best_flag = False
                                    if ms2_match == "":
                                        ms2_match = '\t'
                                best_line = adduct+'\t'+ ms2_match[0:-1]+'\t'+name+ms_score+inchi+formula
                                best_line += "\t" + dmz + score + isotope_match               
                                output.append(i+"\t"+mz+rt+ms2_exp[0:-1]+"\t"+best_line+'\t'+compound+'\n')         
    with open("ms2_output.csv","w")as outf:
        i_d = ""
        for i in output:
            if i.split('\t')[0] == i_d:
                score = float(i.split("\t")[7].replace("","0"))
                if score in best.keys():
                    best[score].append(i)
                else:
                    best[score] = [i]    
            elif i_d == "":
                if i[0] == "I":
                    outf.write(i)
                else:
                    i_d =i.split('\t')[0]
                    best = {}
                    score = float(i.split("\t")[7].replace("","0"))
                    best[score] = [i]    
            else:
                i_d = i.split('\t')[0]
                for b in sorted(best,reverse = True):
                   for o in best[b]:
                       outf.write(o)
                       break
                   break
                best = {}
                score = float(i.split("\t")[7].replace("","0"))
                best[score] = [i]           
    return count

def get_formulas():
    f_c_matches = {}
    try:
        file =  open('compound_dict.txt','r')
        f_c_matches = json.load(file)
        file.close()
    except:
        pass
    output = list()
    files = os.listdir(".\\MAT1")
    for fl in files:
        if fl.endswith('fgt'):
            file = fl.split(".")[0]
            with open(".\\MAT1\\"+file+".fgt",'r') as f_file:
                i = file.split('_')[0] + "\t"
                adduct = file.split('_')[1]
                rt = ""
                mz = ""
                dmz = ""
                score = ""
                isotope_source = ""
                isotope_match = ""
                m1 = ''
                m2 = ''
                compound = '\t'
                with open(".\\MAT1\\"+file + ".MAT", 'r') as mat:
                    for line in mat:
                        if line[:15] == "RETENTIONTIME: ":
                            rt = line[15:-1] + "\t"
                entries = f_file.read().split("\n\n")
                count = 0
                for entry in entries:
                    if count ==1:
                        break
                    else:
                        count += 1
                        for line in entry.split('\n'):
                            if line[:5] == 'NAME:':
                                formula = (line[6:])
                                if formula in f_c_matches.keys():
                                    compound = str(f_c_matches[formula])
                                formula = formula + "\t"
                            if line[:10] == 'EXACTMASS:':
                                mz = line[11:-1]+ "\t"
                            if line[:15] == 'MASSDIFFERENCE:':
                                dmz = line[16:-1]+"\t"
                            if line [:11] == "TOTALSCORE:":
                                score = line[12:-1] + "\t"
                                #if float(score) > 2.5:
                                
                            if line[:24] == "ISOTOPICINTENSITY[M+1]: ":
                                isotope_match = "1;" + line[24:]
                                m1 = float(line[24:])
                                
                            if line[:24] == "ISOTOPICINTENSITY[M+2]: ":
                                isotope_match = isotope_match + ";" + line[24:] + "\t"
                                m2 = float(line[24:])
                            if line[:19] == "ISOTOPICDIFF[M+1]: ":
                                isotope_diff =  float(line[19:])
                                m_d1 = round(m1 - isotope_diff,4)
                                isotope_source = "1;" + str(m_d1)
                            if line[:19] == "ISOTOPICDIFF[M+2]: ":
                                isotope_diff =  float(line[19:])
                                m_d2 = round(m2 - isotope_diff,4)
                                isotope_source =isotope_source+ ";" + str(m_d2) + "\t"
                                output.append(i + mz + rt + formula + dmz + score + isotope_source + isotope_match+adduct+'\t'+compound+"\n")
                        formula = ""
                        mz = ""
                        dmz = ""
                        rt = ""
                        score = ""
    with open('compound_formulas.txt','w') as outf:
        outf.write('index\tmz\trt\tformula\tdmz\tscore\tisotopes_source\tisotopes_match\tAdduct\tCompounds\textra\n')
        i_d = ""
        best = ""
        for i in output:
            if i.split('\t')[0] == i_d:
                i_score = float(i.split("\t")[5].replace("","0"))
                best_score = float(best.split("\t")[5].replace("","0"))
                if i_score > best_score:
                    best = i      
            elif i_d == "":
                if i[0] == "I":
                    outf.write(i)
                else:
                    i_d =i.split('\t')[0]
                    best = i    
            else: 
               outf.write(best)
               i_d = i.split('\t')[0]
               best = i                
    return 1