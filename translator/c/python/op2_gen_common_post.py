import os  
import re
import op2_gen_common

def gen_user_headfile(kernel_names, path):
	#读取路径下的.h和.cpp文件
	datanames = os.listdir(path)
	list = []
	for i in datanames:
		list.append(i)

	list_h=[]
	list_cpp=[]
	list_h_cpp=[]
	for i in list:
		if(i[-2:]=='.h' and i[:-2] not in kernel_names):
			list_h.append(i)
			list_h_cpp.append(i)
		elif(len(i)>4 and i[-4:]=='.cpp'):
			list_cpp.append(i)
			list_h_cpp.append(i)
	return (list_h,list_cpp,list_h_cpp)

def match_funtion_h( path,list_h):
	#找出文件列表中list_h匹配function_return_type的值(函数声明)
	function_return_type = r'''
                    (\s*) #匹配所有的空白字符
                    (   
                       (const)?(volatile)?(static)?\s*(inline)?\s*(extern)?\s*
                                (
                                    (VOID)|(void)|(enum)|
                                        ((unsigned)?(signed)?(long)?\s*(int)|(char)|(float)|(short)|(long)|(double))|
                                        (bool)|(struct\s*\w+)|(union\s*\w+)|(wait_queue_t)|(wait_queue_head_t)
                                )
                                \s+(fastcall)?
                        
                    )  #匹配函数返回类型
                    (\s*(\*)?\s*) #识别有无指针类型*，以及空白字符
                    (\w+) #识别函数名称 
                    ((\s*)(\()) #识别函数开始小括号          
'''
	pat1 = re.compile(function_return_type, re.X)
	funtions_name_dirt={}
	funtions_name_set=set()
	for file_name in list_h:
		#funtions_name.append(file_name)
		path1=path+file_name
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()
		all_the_text = op2_gen_common.comment_remover(all_the_text)
		ret = pat1.finditer(all_the_text) #包含
		for i in ret:
			#funtions_name.append(all_the_text[i.span()[0]:i.span()[1]].split(" ")[-1][:-1].strip())
			funtions_name_dirt[i.groups()[-4]]=file_name
			funtions_name_set.add(i.groups()[-4])
	return funtions_name_dirt,funtions_name_set

def match_funtion_kernel( path,kernel_names,funtions_name_set):
	#匹配kernel中调用的函数
	function_return_type2 = r'''
                    (\s*) #匹配所有的空白字符
                    (\w+)  #匹配函数返回类型
                    ((\s*)(\()) #识别函数开始小括号          
'''
	pat2 = re.compile(function_return_type2, re.X)
	kernel_call_funtion=set()
	for kernel_name in kernel_names:
		path1 = path+kernel_name+'.h'
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()
		all_the_text = op2_gen_common.comment_remover(all_the_text)
		ret = pat2.finditer(all_the_text)
		for i in ret:
			#print(i.group(2))
			if i.group(2) in funtions_name_set:
				kernel_call_funtion.add(i.group(2))
	return kernel_call_funtion

def loop_match_funtion(kernel_call_funtion,path,funtions_name_dirt,funtions_name_set):
	len_kernel_call_funtion=0
	file_text=''
	while len_kernel_call_funtion!=len(kernel_call_funtion):
		len_kernel_call_funtion=len(kernel_call_funtion)
		file_text=''
		for funtion_name in kernel_call_funtion:
			#读取文件
			f = open(path+funtions_name_dirt[funtion_name], 'r',encoding='raw_unicode_escape')
			kernel_text = f.read()
			f.close()
			#去掉注释
			kernel_text = op2_gen_common.comment_remover(kernel_text)
			kernel_text = op2_gen_common.remove_trailing_w_space(kernel_text)
			#提取文件中函数实现
			function_return_type = r'''
					(\s*) #匹配所有的空白字符
					(   
						(
							(VOID)|(void)|(enum)|
							((unsigned)?(signed)?(long)?\s*(int)|(char)|(float)|(short)|(long)|(double))|
							(bool)|(struct\s*\w+)|(union\s*\w+)|(wait_queue_t)|(wait_queue_head_t)
						)
						\s+(fastcall)?
					)  #匹配函数返回类型
					(\s*(\*)?\s*) #识别有无指针类型*，以及空白字符
					('''+funtion_name +''') #识别函数名称 
					((\s*)(\()) #识别函数开始小括号          
		'''
			p = re.compile(function_return_type, re.X)
			i = p.search(kernel_text).start()
			i2 = i
			j = kernel_text[i:].find('{')
			jj = kernel_text[i:].find(';')
			if j==-1 or (jj!=-1 and jj<j):
				f = open(path+funtions_name_dirt[funtion_name][:-2]+'.cpp', 'r',encoding='raw_unicode_escape')
				kernel_text = f.read()
				f.close()
				kernel_text = op2_gen_common.comment_remover(kernel_text)
				kernel_text = op2_gen_common.remove_trailing_w_space(kernel_text)
				function_return_type = r'''
						(\s*) #匹配所有的空白字符
						(   
								(
								(VOID)|(void)|(enum)|
									((unsigned)?(signed)?(long)?\s*(int)|(char)|(float)|(short)|(long)|(double))|
									(bool)|(struct\s*\w+)|(union\s*\w+)|(wait_queue_t)|(wait_queue_head_t)
								)
								\s+(fastcall)?
						)  #匹配函数返回类型
						(\s*(\*)?\s*) #识别有无指针类型*，以及空白字符
						('''+funtion_name +''') #识别函数名称 
						((\s*)(\()) #识别函数开始小括号          
				'''
				p = re.compile(function_return_type, re.X)
				i = p.search(kernel_text).start()
				i2 = i
				j = kernel_text[i:].find('{')
			k = op2_gen_common.para_parse(kernel_text, i+j, '{', '}')
			signature_text = kernel_text[i:i+j]
			l = signature_text[0:].find('(')
			head_text = signature_text[0:l].strip() #save function name
			m = op2_gen_common.para_parse(signature_text, 0, '(', ')')
			signature_text = signature_text[l+1:m]
			body_text = kernel_text[i+j+1:k]

			file_text='__device__' + ' ' + head_text + '(' + signature_text +');\n' + file_text +'__device__' + ' ' + head_text + '(' + signature_text +')\n{' + body_text +'\n}\n'

		#查找kernel中调用的函数
		function_return_type2 = r'''
					(\s*) #匹配所有的空白字符
					(\w+)  #匹配函数返回类型
					((\s*)(\()) #识别函数开始小括号          
		'''
		pat2 = re.compile(function_return_type2, re.X)
		kernel_call_funtion=set()
		ret = pat2.finditer(file_text)
		for i in ret:
			if i.group(2) in funtions_name_set:
				kernel_call_funtion.add(i.group(2))
	return file_text,kernel_call_funtion

def handle_text(path,file_text,kernel_call_funtion,list_cpp,list_h):
	#替换函数名
	for funtion_name in kernel_call_funtion:
		file_text = re.sub('\\b'+funtion_name+'\\b', funtion_name+'_gpu',file_text)

	#注释getchar()
	file_text = re.sub('getchar\(\);', r'/*getchar();*/',file_text)

	# 匹配#define
	function_return_type = r'''
			(\s*) #匹配所有的空白字符
			(\#define)
			(\s+)
			(\w+)
			(\s+)
			(\S+)
	'''
	pat1 = re.compile(function_return_type, re.X)
	define_dirt={}
	for file_name in list_cpp:   
		path1=path+file_name
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()    
		all_the_text = op2_gen_common.comment_remover(all_the_text)
		ret = pat1.finditer(all_the_text) #包含
		for i in ret:
			define_dirt[i.group(4)]=i.group(6)
	for define_name in define_dirt:
		temp=''
		temp+='#ifndef ' + define_name + '\n'
		temp+='#define '+ define_name + ' ' +define_dirt[define_name] + '\n'
		temp+='#endif\n'
		file_text=temp+file_text

	#引入头文件
	for hname in list_h:
		file_text='#include \"../'+ hname + '\"\n' + file_text

	#写入文件
	fid = open('./cuda/user_cuda.h','w')
	fid.write(file_text)
	fid.close()

def replace_kernel(kernel_names,kernel_call_funtion):
	for kernel_name in kernel_names:
		path1='./cuda/'+kernel_name+'_kernel.cu'
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()
		for funtion_name in kernel_call_funtion:
			all_the_text = re.sub('\\b'+funtion_name+'\\b', funtion_name+'_gpu',all_the_text)
		fid = open(path1,'w')
		fid.write(all_the_text)
		fid.close()

def op2_gen_cuda_post(kernel_names):
	path='./'
	list_h,list_cpp,list_h_cpp = gen_user_headfile(kernel_names,path)

	#匹配所有头文件中包含的函数
	funtions_name_dirt,funtions_name_set = match_funtion_h(path,list_h)

	#匹配所有kernel文件中包含的函数
	kernel_call_funtion = match_funtion_kernel( path,kernel_names,funtions_name_set)

	#将找到的函数写入file_text
	file_text,kernel_call_funtion = loop_match_funtion(kernel_call_funtion,path,funtions_name_dirt,funtions_name_set)

	#写入头文件
	handle_text(path,file_text,kernel_call_funtion,list_cpp,list_h)

	replace_kernel(kernel_names,kernel_call_funtion)


def op2_gen_openmp_post(kernel_names):
	path='./'
	list_h,list_cpp,list_h_cpp = gen_user_headfile(kernel_names,path)
	file_text=''
	for hname in list_h:
		file_text='#include \"../'+ hname + '\"\n' + file_text
	#写入文件
	fid = open('./openmp/user_openmp.h','w')
	fid.write(file_text)
	fid.close()

def op2_gen_openmp_extern(kernel_names , ex_name):
	path='./'
	list_h,list_cpp,list_h_cpp = gen_user_headfile(kernel_names,path)
	global_pattern = r'const[\s]+[\w]+[\s]+'+ex_name
	
	for file_name in list_h_cpp:
		path1=path+file_name
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()    
		all_the_text = op2_gen_common.comment_remover(all_the_text)
		if re.findall(global_pattern, all_the_text):
			return True
	return False

def op2_gen_genseq_post(kernel_names):
	path='./'
	list_h,list_cpp,list_h_cpp = gen_user_headfile(kernel_names,path)
	file_text=''
	for hname in list_h:
		file_text='#include \"../'+ hname + '\"\n' + file_text
	#写入文件
	fid = open('./seq/user_genseq.h','w')
	fid.write(file_text)
	fid.close()
	
struct_def_dirt={}
def gen_user_datatype(kernel_names):
	path='./'
	list_h,list_cpp,list_h_cpp = gen_user_headfile(kernel_names,path)

	function_return_type = r'''
			(\s*) #匹配所有的空白字符
			(typedef)  #匹配函数返回类型
			(\s*) #匹配所有的空白字符
			(struct)
	'''
	pat1 = re.compile(function_return_type, re.X)
	struct_def=set()
	for file_name in list_h:
		path1=path+file_name
		file_object = open(path1,'r',encoding='raw_unicode_escape')
		all_the_text = file_object.read()
		file_object.close()
		all_the_text = op2_gen_common.comment_remover(all_the_text)
		ret = pat1.finditer(all_the_text)
		#print(all_the_text)
		for match_i in ret:
			i=match_i.start()
			j = all_the_text[i:].find('{')
			if j==-1:
				print('find struct error in ',file_name,':',i)
				exit(1)
			k = op2_gen_common.para_parse(all_the_text, i+j, '{', '}')
			jj = all_the_text[k:].find(';')
			struct_def.add(all_the_text[k+1:k+jj].strip())
			struct_def_dirt[all_the_text[k+1:k+jj].strip()]=all_the_text[i+j:k+1]
	print(kernel_names)
	file_text=''
	for i in struct_def:
		file_text+='inline int type_error(const '+ i +' *a, const char *type)\n'
		file_text+='{\n'
		file_text+='\t(void)a;\n'
		file_text+='\treturn (strcmp(type, "'+ i +'") && strcmp(type, "'+ i +':soa"));\n'
		file_text+='}\n\n'
	file_text+='template <class T>\n'
	file_text+='inline int type_error(const T *a, const char *type)\n'
	file_text+='{\n'
	file_text+='\treturn type_error(*a, type);\n'
	file_text+='}\n\n'

	op2_path=os.environ.get('OP2_INSTALL_PATH')
	path2=op2_path+'/c/include/op_user_datatypes.h'
	fid = open(path2,'w')
	fid.write(file_text)
	fid.close()

def gen_user_datatype_mpi(macro_defs):
	file_text=''' 
	template <typename T>
	void gather_data_hdf5(op_dat dat, char *usr_ptr, int low, int high);

	template <typename T, void (*F)(FILE *, int, int, T *, const char *)>
	void write_file(op_dat dat, const char *file_name);

	template <typename T, const char *fmt>
	void write_txt(FILE *fp, int g_size, int elem_size, T *g_array,const char *file_name);

	template <typename T>
	void write_bin(FILE *fp, int g_size, int elem_size, T *g_array,const char *file_name);

	'''
	for i in macro_defs:
		file_text+='\n#ifndef '+i+'\n'
		file_text+='#define '+i+' '+macro_defs[i]+'\n'
		file_text+='#endif'+'\n'
	
	for i in struct_def_dirt:
		file_text += '\ntypedef struct'+struct_def_dirt[i]+i+';\n'
		file_text += 'MPI_Datatype MPI_'+i.upper()+';\n'


	struct_re=''
	for struct_name in struct_def_dirt:
		struct_re+='|('+struct_name+')'

	not_used_struct=[]
	for struct_name in struct_def_dirt:
	#     print(struct_name)
		function_return_type = r'''
				(\s*) #匹配所有的空白字符
				(unsigned)? #匹配函数返回类型
				(\s*) #匹配所有的空白字符
				((int)|(char)|(float)|(double)'''+struct_re+''')
				(\s*) #匹配所有的空白字符
				(\*)?
				(\w+)
				(.*)
				(;)
		'''
		pat1 = re.compile(function_return_type, re.X)
		all_the_text=struct_def_dirt[struct_name]
		ret = pat1.finditer(all_the_text)
		
		struct_i=[]
		
		for i in ret:
	#         print(i.groups()[3],i.groups()[-4],i.groups()[-3],i.groups()[-2])
			struct_i.append([i.groups()[3],i.groups()[-4],i.groups()[-3],i.groups()[-2]])
	#         print(i.groups())
		if len(struct_i) == 0:
			not_used_struct.append(struct_name)
			continue
		#确定数据类型
		struct_i_lenght=str(len(struct_i))
		struct_i_text= '\nvoid struct_' + struct_name +'()\n{\n'
		struct_i_text += 'MPI_Datatype array_of_types['+struct_i_lenght+'] = {'
		for i in struct_i:
			if i[1]=='*':
				struct_i_text += 'MPI_DOUBLE,'
			else:
				struct_i_text += 'MPI_'+i[0].upper()+',';
		struct_i_text = struct_i_text[:-1]+'};\n'
		
		#确定偏移量
		struct_i_text += struct_name + ' ' + struct_name + '_user;\n'
		struct_i_text += 'MPI_Aint array_of_displacements[' + struct_i_lenght +'];\n'
		for i in range(len(struct_i)):
			struct_i_text += 'MPI_Get_address(&' + struct_name + '_user' + '.' + struct_i[i][2] +', &array_of_displacements[' + str(i) + ']);\n'
		
		
		for i in range(len(struct_i),0,-1):
			struct_i_text += 'array_of_displacements[' + str(i) + '] = ' +'array_of_displacements[' + str(i) + '] - array_of_displacements[0];\n'
		struct_i_text += 'array_of_displacements[0] = 0;\n'

		#确定数组大小
		struct_i_text += 'int array_of_blocklengths[' + struct_i_lenght + '] = { '
		for i in struct_i:
			if i[3]=='':
				struct_i_text += str(1)+', '
			else:
				macro_to_int=''
				for c in i[3]:
					macro_to_int = macro_to_int.strip()
					if c==']':
						if macro_to_int.isdigit():
							struct_i_text += macro_to_int+'*'
						elif macro_to_int in macro_defs:
							struct_i_text += macro_defs[macro_to_int]+'*'
						else:
							print(struct_name + ' macro error!')
							exit()
					elif c=='[':
						macro_to_int=''
					else:
						macro_to_int+=c
				struct_i_text = struct_i_text[:-1]+', '
		struct_i_text = struct_i_text[:-2]+' };\n'
		
		#创建MPI结构体
		struct_i_text += 'MPI_Type_create_struct(' + struct_i_lenght + ', array_of_blocklengths, array_of_displacements, array_of_types, &'
		struct_i_text += 'MPI_' + struct_name.upper() +');\n'
		struct_i_text += 'MPI_Type_commit(&' + 'MPI_' + struct_name.upper() +');\n}\n'
		
		#重载_mpi_allgather，_mpi_gather
		struct_i_text += '\nvoid _mpi_allgather('+ struct_name + '*l, '+ struct_name +' *g, int size, int *recevcnts, int *displs, MPI_Comm comm)\n'
		struct_i_text += '{\n'
		struct_i_text += 'struct_' + struct_name +'();\n'
		struct_i_text += 'MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_' + struct_name.upper() + ', comm);\n}\n'
		
		struct_i_text += '\nvoid _mpi_gather('+ struct_name +' *l, '+ struct_name +' *g, int size, int *recevcnts, int *displs, MPI_Comm comm)\n'
		struct_i_text += '{\n'
		struct_i_text += 'struct_' + struct_name +'();\n'
		struct_i_text += 'MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_' + struct_name.upper() + ',MPI_ROOT, comm);\n}\n'
		
		file_text += struct_i_text

	file_text += '\nvoid fetch_data_hdf5_user_datatype(op_dat dat, char *usr_ptr, int low, int high)\n'
	file_text += '{\n'

	for struct_name in struct_def_dirt:
		if struct_name not in not_used_struct:
			file_text += 'if (strcmp(dat->type,\"' + struct_name + '\") == 0)\n'
			file_text += 'gather_data_hdf5<' + struct_name + '>(dat, usr_ptr, low, high);\nelse '

	file_text += '\nprintf(\"Unknown type %s, cannot error in fetch_data_hdf5() \\n\", dat->type);\n}\n'

	file_text += '\nextern const char fmt_struct[] = \"%f \";\n\n'

	file_text += '\nvoid print_dat_to_txtfile_mpi_user_datatypes(op_dat dat, const char *file_name)\n{\n'

	for struct_name in struct_def_dirt:
		if struct_name not in not_used_struct:
			file_text += 'if (strcmp(dat->type,\"' + struct_name + '\") == 0)\n'
			file_text += 'write_file<' + struct_name + ', write_txt<' + struct_name + ', fmt_struct>>(dat, file_name);\nelse '

	file_text += '\nprintf(\"Unknown type %s, cannot be written to file %s\\n\", dat->type,file_name);\n}\n\n'

	file_text += '\nvoid print_dat_to_binfile_mpi_user_datatypes(op_dat dat, const char *file_name)\n{\n'

	for struct_name in struct_def_dirt:
		if struct_name not in not_used_struct:
			file_text += 'if (strcmp(dat->type,\"' + struct_name + '\") == 0)\n'
			file_text += 'write_file<' + struct_name + ', write_bin<' + struct_name + '>>(dat, file_name);\nelse '

	file_text += '\nprintf(\"Unknown type %s, cannot be written to file %s\\n\", dat->type,file_name);\n}\n\n'

	op2_path=os.environ.get('OP2_INSTALL_PATH')
	path2=op2_path + '/c/src/mpi/op_mpi_user_datatypes.h'
	fid = open(path2,'w')
	fid.write(file_text)
	fid.close()