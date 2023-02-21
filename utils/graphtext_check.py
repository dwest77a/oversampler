



def import_graph(ecv_type, label, daynight, years, months, family):
    data = []  
    base_path = '/home/users/dwest77/Documents/ECV_Images/output_files/graph_texts/{}'.format(ecv_type+label)
    output_file = '/{}_{}_{}_{}_{}_{}g.txt'.format(family, ecv_type+label, daynight, 'sea', years, months)
        
    f = open(base_path+output_file, 'r')
    outstring = f.readlines()
    header = outstring[0].replace(' ','').split(',')
    for index in range(1,len(outstring)):
        line = outstring[index].replace('\n','')
        line = line.replace(' ','').split(',')
           
        data.append(float(line[1]))
       
        
        f.close()	
    return data
    
    
cris_nh3_day = import_graph('nh3','','day','20162018','112','cris')
cris_nh3_night = import_graph('nh3','','night','20162018','112','cris')
#iasi_nh3_day = import_graph('nh3','','day','20082018','112','iasi')
#iasi_nh3_night = import_graph('nh3','','night','20162018','112','iasi')

for index in range(len(cris_nh3_day)):
    print(cris_nh3_day[index] - cris_nh3_night[index])
    
