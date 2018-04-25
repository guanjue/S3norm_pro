data = open('raw_5p_rc_list.txt','r')
data_dict = {}

for records in data:
	tmp_info = records.split()
	ct = tmp_info[0].split('.')[0]
	if ct in data_dict:
		data_dict[ct].append(tmp_info[0])
	else:
		data_dict[ct] = [tmp_info[0]]

data.close()

result = open('raw_5p_rc_list_2d.txt','w')

for ct in data_dict:
	l = data_dict[ct]
	if (len(l)==2):
		result.write(l[0] + '\t' + l[1] + '\n')
	elif (len(l)==3):
		result.write(l[0] + '\t' + l[1] + '\n')
		result.write(l[1] + '\t' + l[2] + '\n')

result.close()

