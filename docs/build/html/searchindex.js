Search.setIndex({docnames:["gaddlemaps","gaddlemaps.components","gaddlemaps.parsers","index","modules","overview"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.index":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["gaddlemaps.rst","gaddlemaps.components.rst","gaddlemaps.parsers.rst","index.rst","modules.rst","overview.rst"],objects:{"":{gaddlemaps:[0,0,0,"-"]},"gaddlemaps.Alignment":{SIGMA_SCALE:[0,2,1,""],STEPS_FACTOR:[0,2,1,""],align_molecules:[0,3,1,""],end:[0,3,1,""],exchange_map:[0,2,1,""],init_exchange_map:[0,3,1,""],start:[0,3,1,""],write_comparative_gro:[0,3,1,""]},"gaddlemaps.Chi2Calculator":{chi2_molecules:[0,3,1,""]},"gaddlemaps.ExchangeMap":{equivalences:[0,3,1,""]},"gaddlemaps.Manager":{add_end_molecule:[0,3,1,""],add_end_molecules:[0,3,1,""],align_molecules:[0,3,1,""],calculate_exchange_maps:[0,3,1,""],complete_correspondence:[0,3,1,""],extrapolate_system:[0,3,1,""],from_files:[0,3,1,""],molecule_correspondence:[5,2,1,""],parse_restrictions:[0,3,1,""]},"gaddlemaps.components":{Atom:[1,1,1,""],AtomGro:[1,1,1,""],AtomTop:[1,1,1,""],Molecule:[1,1,1,""],MoleculeTop:[1,1,1,""],Residue:[1,1,1,""],System:[1,1,1,""],SystemGro:[1,1,1,""],are_connected:[1,4,1,""]},"gaddlemaps.components.Atom":{atom_gro:[1,3,1,""],atom_top:[1,3,1,""],copy:[1,3,1,""],gro_resid:[1,3,1,""],top_resid:[1,3,1,""]},"gaddlemaps.components.AtomGro":{atomid:[1,2,1,""],copy:[1,3,1,""],element:[1,3,1,""],gro_line:[1,3,1,""],name:[1,2,1,""],position:[1,2,1,""],resid:[1,2,1,""],residname:[1,3,1,""],resname:[1,2,1,""],velocity:[1,2,1,""]},"gaddlemaps.components.AtomTop":{bonds:[1,2,1,""],closest_atoms:[1,3,1,""],connect:[1,3,1,""],copy:[1,3,1,""],residname:[1,3,1,""]},"gaddlemaps.components.Molecule":{atoms:[1,3,1,""],atoms_ids:[1,3,1,""],atoms_positions:[1,3,1,""],atoms_velocities:[1,3,1,""],bonds_distance:[1,3,1,""],copy:[1,3,1,""],deep_copy:[1,3,1,""],distance_to:[1,3,1,""],distance_to_zero:[1,3,1,""],from_files:[1,3,1,""],geometric_center:[1,3,1,""],index:[1,3,1,""],molecule_top:[1,3,1,""],move:[1,3,1,""],move_to:[1,3,1,""],remove_atom:[1,3,1,""],resid:[1,3,1,""],residname:[1,3,1,""],resids:[1,3,1,""],residues:[1,3,1,""],resname:[1,3,1,""],resnames:[1,3,1,""],rotate:[1,3,1,""],update_from_molecule_top:[1,3,1,""],write_gro:[1,3,1,""],x:[1,3,1,""],y:[1,3,1,""],z:[1,3,1,""]},"gaddlemaps.components.MoleculeTop":{copy:[1,3,1,""],index:[1,3,1,""],resids:[1,3,1,""],resname_len_list:[1,3,1,""],resnames:[1,3,1,""]},"gaddlemaps.components.Residue":{atoms:[1,3,1,""],atoms_ids:[1,3,1,""],atoms_positions:[1,3,1,""],atoms_velocities:[1,3,1,""],copy:[1,3,1,""],distance_to:[1,3,1,""],distance_to_zero:[1,3,1,""],geometric_center:[1,3,1,""],move:[1,3,1,""],move_to:[1,3,1,""],remove_atom:[1,3,1,""],resid:[1,3,1,""],residname:[1,3,1,""],resname:[1,3,1,""],rotate:[1,3,1,""],update_from_molecule_top:[1,3,1,""],write_gro:[1,3,1,""],x:[1,3,1,""],y:[1,3,1,""],z:[1,3,1,""]},"gaddlemaps.components.System":{add_ftop:[1,3,1,""],add_molecule_top:[1,3,1,""],composition:[1,3,1,""],fgro:[1,3,1,""]},"gaddlemaps.components.SystemGro":{box_matrix:[1,3,1,""],comment_line:[1,3,1,""],composition:[1,3,1,""],molecules_info_ordered_all:[1,3,1,""],molecules_resname_len_index:[1,3,1,""],n_atoms:[1,3,1,""]},"gaddlemaps.parsers":{CoordinatesParser:[2,1,1,""],GroFile:[2,1,1,""],ParserManager:[2,1,1,""],ParserRegistered:[2,1,1,""],dump_lattice_gro:[2,4,1,""],extract_lattice_gro:[2,4,1,""],open_coordinate_file:[2,4,1,""]},"gaddlemaps.parsers.CoordinatesParser":{EXTENSIONS:[2,2,1,""],box_matrix:[2,3,1,""],close:[2,3,1,""],comment:[2,3,1,""],natoms:[2,3,1,""],next:[2,3,1,""],seek_atom:[2,3,1,""],writeline:[2,3,1,""],writelines:[2,3,1,""]},"gaddlemaps.parsers.GroFile":{COORD_START:[2,2,1,""],DEFAULT_COMMENT:[2,2,1,""],DEFAULT_POSTION_FORMAT:[2,2,1,""],EXTENSIONS:[2,2,1,""],NUMBER_FIGURES:[2,2,1,""],box_matrix:[2,3,1,""],close:[2,3,1,""],comment:[2,3,1,""],determine_format:[2,3,1,""],name:[2,3,1,""],natoms:[2,3,1,""],next:[2,3,1,""],parse_atomline:[2,3,1,""],parse_atomlist:[2,3,1,""],position_format:[2,3,1,""],readline:[2,3,1,""],readlines:[2,3,1,""],seek_atom:[2,3,1,""],validate_string:[2,3,1,""],writeline:[2,3,1,""],writelines:[2,3,1,""]},"gaddlemaps.parsers.ParserManager":{parsers:[2,2,1,""],register:[2,3,1,""]},"gaddlemaps.parsers.ParserRegistered":{mro:[2,3,1,""],register:[2,3,1,""]},gaddlemaps:{Alignment:[0,1,1,""],Chi2Calculator:[0,1,1,""],ExchangeMap:[0,1,1,""],Manager:[0,1,1,""],accept_metropolis:[0,4,1,""],calcule_base:[0,4,1,""],comparate_alignment:[0,4,1,""],components:[1,0,0,"-"],find_atom_random_displ:[0,4,1,""],guess_protein_restrains:[0,4,1,""],guess_residue_restrains:[0,4,1,""],interactive_restrictions:[0,4,1,""],minimize_molecules:[0,4,1,""],move_mol_atom:[0,4,1,""],parsers:[2,0,0,"-"],remove_hydrogens:[0,4,1,""],rotation_matrix:[0,4,1,""]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","attribute","Python attribute"],"3":["py","method","Python method"],"4":["py","function","Python function"]},objtypes:{"0":"py:module","1":"py:class","2":"py:attribute","3":"py:method","4":"py:function"},terms:{"3x3":[1,2,5],"abstract":2,"case":[0,1],"class":[0,1,2,3],"default":[0,1,2,5],"final":[0,2,5],"float":[0,1,2,5],"function":[0,1,2,5],"import":[0,5],"int":[0,1,2,5],"new":[0,1,2,5],"return":[0,1,2,5],"static":[0,2],"true":[0,1,2,5],"try":0,ARE:[0,5],For:[0,1,5],Not:[0,5],One:0,THE:[0,5],THe:2,The:[0,1,2,5],Then:[0,5],These:[0,1],With:[1,5],_compar:0,_compon:[0,1,5],_components_top:[1,5],_manag:0,_residu:[0,1,5],abc:2,abcmeta:2,about:[1,2],accept:0,accept_metropoli:0,access:[0,1,5],accord:0,accordion:0,account:0,act:1,actual:[1,5],add:[0,1,5],add_end_molecul:[0,5],add_ftop:[1,5],add_molecule_top:[1,5],advic:[1,5],affect:[1,5],afford:1,after:0,again:0,algorithm:0,align:[0,1,5],align_molecul:[0,5],alignemnt:[0,5],aling_molecul:0,alist:2,all:[0,1,2,5],allow:[0,2,5],almost:0,alphabet:1,alreadi:2,also:[0,1,2,5],although:[0,1,5],amplitud:0,analyz:2,anchura:0,angl:0,ani:[1,5],app_point:0,appear:1,appli:[0,1,5],applic:0,appropri:[1,2],are_connect:1,argument:0,around:[0,1,5],arrai:[0,1,2,5],articl:0,assign:0,associ:[0,5],atom:[0,1,2,5],atom_gro:1,atom_index:0,atom_po:0,atom_top:1,atomgro:[1,5],atomid:1,atomist:0,atomitp:1,atomlin:2,atomlist:2,atoms_id:[1,5],atoms_info:1,atoms_po:0,atoms_posit:[1,5],atoms_veloc:[1,5],atomtop:1,attach:[0,5],attr:2,attribut:[0,1,2,5],auto_guess_protein_restrict:0,autodetect:2,automat:[0,2],avail:1,avoid:[0,1,5],axi:0,back:0,base:[0,1,2],basi:0,basic:[0,1,3],befor:[0,2],begin:2,behav:1,behaviour:2,belong:0,better:[0,1,5],between:[0,1,5],bmimaa_new:0,bmimaa_target:0,bmimcg_new:0,bmimcg_ref:0,bond:[0,1,5],bonded_atom:1,bonds_dist:[1,5],bonds_info:0,bool:[0,1,2,5],both:[0,1,5],boundari:[1,5],box:[0,1,2,5],box_matrix:[1,2],box_vect:[1,5],build:[0,5],built:[0,5],calcul:[0,5],calculate_exchange_map:[0,5],calcule_bas:0,call:[0,2,5],can:[0,1,2,5],carlo:0,center:[0,1,5],certain:0,chang:[1,3,4,5],charact:[1,2],check:[0,1,5],chi2:0,chi2_molecul:0,chi2calcul:0,chosen:[0,5],classmethod:[0,1,2,5],clockwis:0,close:[0,2,5],closest:0,closest_atom:1,closest_r2_atoms_index:0,cls:2,coars:0,code:[1,5],collaps:0,collect:[0,1,5],collinear:0,combin:[0,1,5],come:2,comment:[1,2],comment_lin:1,common:2,compar:0,comparate_align:0,compat:[1,2,5],complet:[0,1,5],complete_correspond:[0,5],complex:[0,1,5],compon:[0,3,4,5],composit:[1,5],compress:[0,5],comput:[0,1,5],condit:[1,5],configur:[0,5],conform:[1,5],connect:1,consecut:1,consid:[0,1,5],constant:0,constitut:[1,5],constraint:0,contain:[0,1,2],contempl:2,content:[3,4],convert:[0,2],coord_start:2,coordin:[0,1,2,5],coordinatespars:2,copi:[1,5],correct:[0,1,2,5],correspond:[0,1,2,5],counter:[0,1,5],crash:0,creat:[0,1,2,5],creation:0,current:[0,1],cursor:2,cut:2,data:[0,1,2,5],decim:2,decor:2,deep:[1,5],deep_copi:[1,5],default_com:2,default_postion_format:2,defaultdict:0,defin:[0,1],deform:0,deformation_typ:[0,5],depend:0,descript:0,desir:[0,5],determin:0,determine_format:2,dict:[0,1,2,5],dict_format:2,dictionari:[0,1,2,5],differ:[0,1,5],different_molecul:1,direct:0,directli:2,discret:0,displ:0,displac:[0,1,2,5],displacement_modul:0,distanc:[0,1,5],distance_to:[1,5],distance_to_zero:[1,5],distribut:0,document:0,doe:[0,1,5],done:[0,1,5],dump_lattice_gro:2,dure:0,dynam:[0,1],each:[0,1,2,5],ear:1,eas:0,element:1,els:[0,1],empti:[0,1,2],end:[0,5],energi:0,energy_0:0,energy_1:0,engin:[0,5],enumer:1,equal:2,equival:0,euclidean:[1,5],evalu:0,everi:[0,5],exampl:[0,5],exchang:[0,5],exchange_map:0,exchangemap:0,exclud:1,extens:[1,2],extract:[1,2],extract_lattice_gro:2,extrapol:[0,1,5],extrapolate_system:[0,5],f_system_gro:[0,5],factor:[0,5],fail:0,fals:[0,1,2,5],featur:[1,5],fgro2:0,fgro:[0,1,5],fgro_out:[0,5],figur:2,file:[0,1,5],file_format:1,filegro:1,filenam:2,fill:0,find:[0,1,5],find_atom_random_displ:0,first:[0,1,2,5],fitp:0,fix:0,fname:0,follow:2,forcefield:0,form:[0,1,2,5],format:[0,1,2,5],format_dict:2,forwardref:[1,2,5],found:[0,2,5],fout:[1,5],from:[0,1,2,5],from_fil:[0,1,5],ftop:[0,1,5],fulli:0,functor:0,futur:[0,1,5],gaddl:[0,1],gaddlemap:[3,5],gener:[0,1,5],generet:2,geometr:[0,1,5],geometric_cent:[1,5],given:[0,1,5],global_index:2,goe:0,grain:0,gro:[0,1,2,5],gro_lin:1,gro_resid:1,groaa:0,grofil:[1,2],gromac:[0,1,2,5],guess:[0,2,5],guess_protein:[0,5],guess_protein_restrain:[0,5],guess_residue_restrain:0,hand:[0,1],handi:0,has:[0,1,2,5],hash:1,have:[0,1,5],header:2,help:[1,5],how:[0,1,5],howev:0,hydrogen:[0,5],id_num:1,identifi:[1,5],ids:[1,5],ignore_hydrogen:[0,5],implement:[0,2,5],includ:[0,1,5],independ:[1,5],index:[0,1,2,3,5],individu:[0,5],info:2,inform:[0,1,2,5],inherit:[1,5],init:[0,5],init_exchange_map:0,initi:[0,1,2,5],input:[0,1,2,5],insert:2,instanc:[0,1,5],instead:0,integ:[0,1,2],interactive_restrict:0,interest:0,interfac:0,inv:[1,5],invers:[1,5],ioerror:[0,1,2,5],ipywidget:0,iter:1,itp:[0,1,2,5],itp_top:[],itpaa:0,its:[0,1,2,5],joint:0,jupyt:0,just:[0,1,2,5],kei:[0,1,2,5],kept:2,keyerror:[0,5],know:0,larg:1,larger:0,last:0,lattic:[1,2],lead:2,length:[0,1,5],like:[1,2],line:[0,1,2],linearli:0,link:[1,5],list:[0,1,2,5],list_atomlist:2,load:[0,1,5],local:0,look:[0,5],loop:[0,5],low:0,lower:1,mai:[0,1,5],maintain:1,make:[1,5],manag:[0,1,5],manger:[],mani:[0,1,5],manipul:1,manual:[0,5],map:[0,1],mass:[0,1,5],match:[0,1,5],matrix:[0,1,2,5],mean:[1,5],memori:1,method:[0,1,2,5],minim:[0,5],minimize_molecul:0,miss:[1,5],mobil:0,mode:2,modifi:[0,1,2,5],modified_atoms_po:0,modul:[3,4,5],mol1:0,mol1_posit:0,mol2:0,mol2_bonds_info:0,mol2_com:0,mol2_posit:0,mol_index:2,mol_top:[1,5],molec:0,molecul:[1,3,4,5],molecular:[0,2,5],molecule_correspond:[0,5],molecule_top:[1,5],moleculegro:0,molecules_info_ordered_al:1,molecules_resname_len_index:1,moleculetop:[1,5],mont:0,more:[0,1,2,5],most:[0,1,5],move:[0,1,5],move_mol_atom:0,move_to:[1,5],mro:2,mtop:[1,5],multipl:[0,5],must:[0,1,5],n_atom:1,n_step:0,name:[0,1,2,5],natom:[1,2],ndarrai:[0,1,2,5],nearest:0,necessari:0,need:1,new_atom:1,new_molecule_gro:[1,5],new_posit:[0,1,5],new_residu:[1,5],new_restrict:[0,5],new_str:2,next:2,non:1,none:[0,1,2,5],nonetyp:[0,2,5],normal:0,note:[0,1,5],notebook:0,now:[0,1,2,5],number:[0,1,2,5],number_figur:2,number_of_atoms_with_resname_1:1,number_of_atoms_with_resname_2:1,numpi:[0,1,2,5],object:[0,1,2,5],obtain:[0,1],offset1:0,offset2:0,offset:0,old:[1,5],onc:[0,5],one:[0,1,2,5],ones:0,onli:[0,1,2,5],open:[0,2],open_coordinate_fil:2,oper:0,optim:[0,5],optimum:[0,5],option:[0,1,2,5],order:[1,2],origin:[1,5],orthonorm:0,other:[0,1,5],otherwis:0,output:[0,1],over:[0,5],overlap:[0,5],overwrit:[0,5],packag:[3,4,5],page:3,pair:[0,5],parallel:[0,5],paramet:[0,1,2,5],pars:[0,1,2,5],parse_atomlin:2,parse_atomlist:2,parse_restrict:[0,5],parsed_gro_lin:1,parser:[0,1,3,4],parsermanag:2,parserregist:2,part:[1,5],pass:[0,1,5],path:[0,1,2,5],pbc:2,per:0,perform:[0,1,5],period:[1,5],perpendicular:0,pick:0,place:[1,2,5],plane:0,point:[0,1,5],popc:[0,5],pos:0,posit:[0,1,2,5],position_format:2,possibl:[0,2,5],postion:2,prevent:0,previous:[0,5],privileg:[0,5],process:[0,1,5],program:[0,2],properli:2,properti:[0,1,2,5],protein:[0,5],provid:[0,2],publish:0,python:[2,3,4],quit:0,r1_atom_index:0,radian:0,radiu:0,rais:[0,1,2,5],random:0,randomli:0,read:2,read_lin:1,readi:0,readlin:2,refer:[0,1],referenc:[0,5],refmolecul:0,refresh:0,regist:2,regul:0,reimplement:[1,5],relev:0,remain:[0,1,2,5],remov:[0,1,5],remove_atom:[1,5],remove_hydrogen:0,repeat:0,repetit:1,replac:[0,1,5],repres:0,represent:0,request:1,res1:0,res2:0,resid:[1,5],residnam:[1,5],residu:[0,1,5],resnam:[1,2,5],resname_1:1,resname_2:1,resname_len_list:1,resolut:[0,2,5],resolv:[0,5],respect:0,rest:[0,1,5],restrain:[0,5],restrict:[0,5],restriction_widget:0,result:0,right:2,root:[1,5],rotat:[0,1,5],rotation_matrix:[0,1,5],row:[0,2],run:[0,5],safe:1,sai:[0,5],same:[0,1,5],same_com:0,save:[0,5],scale:0,scale_factor:[0,5],search:3,second:[0,1,2,5],see:0,seek_atom:2,self:[0,1,2,5],serv:0,set:[0,1,5],sever:2,shape:1,should:[0,1,2,5],sigma:0,sigma_scal:0,sim_typ:0,simplifi:3,simul:[3,4,5],singl:0,smaller:[0,2],solvent:[0,5],some:[0,2],sourc:[0,1,2,5],space:0,speci:0,specifi:[0,1,5],squar:0,stack:0,start:[0,1,5],step:[0,1,5],steps_factor:0,still:0,storag:1,store:[1,5],str:[0,1,2,5],string:[0,1,2,5],structur:[0,1,5],style:0,subclass:2,submodul:[1,2],subpackag:[3,4],sum:0,summar:1,support:1,sure:[1,5],system:[0,1,2,5],system_compon:[0,5],systemerror:[0,5],systemgro:1,tab:0,take:[0,2],taken:[0,1],targetmolecul:0,than:[0,1,2,5],thei:[0,1,5],them:[0,1,5],theta:0,thi:[0,1,2,5],those:[0,5],three:0,through:[0,1,5],time:[0,1,5],todo:0,toe:[0,5],tool:[0,2],top:[1,5],top_resid:1,topolog:[0,1,5],total:2,transform:0,translat:[0,5],tupl:[0,1,2,5],two:[0,1,5],type:[0,1,2,5],typeerror:[0,1,5],undesir:[1,5],unfavor:0,union:[0,1,2,5],uniqu:1,unitari:0,updat:[1,5],update_from_molecule_top:[1,5],usag:2,use:[1,5],used:[0,1,2,5],useful:[1,2,5],using:[0,5],usual:[0,5],valid:[0,2,5],validate_str:2,valu:[0,1,5],valueerror:[0,1,5],vec1:0,vector:[0,1,2,5],veloc:[1,2,5],veri:[0,1,5],verifi:2,vertic:0,vet_molecul:[0,5],view:0,virtual:2,visual:0,visualiza:0,vte:[0,5],vte_gro:[0,5],vte_itp:[0,5],vte_molecul:[0,5],wai:[0,2],want:[0,1,5],warrant:2,were:[0,2,5],wether:2,when:[0,1,5],where:[0,2,5],whether:2,which:[0,1,5],whole:0,widget:0,without:[0,1,2,5],work:[0,1],workaround:1,worri:1,wrap:[1,5],write:[0,1,2,5],write_comparative_gro:0,write_gro:[1,5],writelin:2,written:2,wrong:[0,1,5],you:[0,1,5],your:[0,5]},titles:["gaddlemaps package","gaddlemaps.components package","gaddlemaps.parsers package","Welcome to Gaddle Maps\u2019s documentation!","gaddlemaps","GADDLE Maps simplified documentation"],titleterms:{"class":5,basic:5,chang:0,compon:1,content:[0,1,2],document:[3,5],file:2,gaddl:[3,5],gaddlemap:[0,1,2,4],indic:3,map:[3,5],modul:[0,1,2],molecul:0,molecular:1,packag:[0,1,2],parser:2,python:0,simplifi:5,simul:[0,1,2],subpackag:0,tabl:3,welcom:3}})