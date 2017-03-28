# Code adapted from Sheng (done on framingham computer)

setwd("~/Framingham Files/work_Joelle")

is_subset = function(s1,s2) {  	
	setequal(s1, intersect(s1,s2))	
}

####check whether the sets of shareids zuoheng and xiaowei retained
####are subsets of the set of shareids in the file "shareid_all"

## All individuals ID 
data_shareid_all=read.table(file="../work_Sheng/shareid_all", header=TRUE, sep="\t", fill=TRUE) #14531 rows
shareid_all = data_shareid_all$shareid
shareid0_all = shareid_all[which(data_shareid_all$IDTYPE==0)] # IDTYPE refers to the cohort
shareid1_all = shareid_all[which(data_shareid_all$IDTYPE==1)]
shareid2_all = shareid_all[which(data_shareid_all$IDTYPE==2)] # Only 102 individuals in cohort 2
shareid3_all = shareid_all[which(data_shareid_all$IDTYPE==3)] 
table(data_shareid_all$IDTYPE)

### Other data set containin related individuals
data_pheno_zuoheng=read.table(file="../work_Sheng/pheno_zuoheng",header=TRUE,sep="\t")  #7678 rows
shareid_related_zhw = data_pheno_zuoheng$shareid
# Break by cohort membership
shareid0_related_zhw = intersect(shareid_related_zhw,shareid0_all)
shareid1_related_zhw = intersect(shareid_related_zhw,shareid1_all)
shareid2_related_zhw = intersect(shareid_related_zhw,shareid2_all)
shareid3_related_zhw = intersect(shareid_related_zhw,shareid3_all)

length(shareid_related_zhw)
length(shareid0_related_zhw)+length(shareid1_related_zhw)+length(shareid2_related_zhw)+length(shareid3_related_zhw) # They do match


##### the set of genotyped individuals
data_indgeno = read.table(file="../work_Sheng/ind_geno",header=TRUE, sep="\t")
shareid_indgeno = data_indgeno$shareid

## Same set of individuals in the two files
shareid_xiaowei_all <- read.table(file="../work_Sheng/ShareID_chr6.txt",header=FALSE)$V1
setequal(shareid_indgeno, shareid_xiaowei_all) # Same set of genotyped individuals

shareid_failxw = read.table(file="../work_Sheng/SidFailQC.txt",header=FALSE,sep="\t")$V1 # People who failed quality control checks
is_subset(shareid_failxw,shareid_xiaowei_all)

## save set of GT individuals who did not fail QC
shareid_xww = setdiff(shareid_xiaowei_all,shareid_failxw) ### contains 8597 individuals

## extract the set of unrelated individuals from xww set
shareid_unrelated_all = shareid_all[which(is.na(data_shareid_all$pedno))]  # 1519 unrelated individuals
shareid_unrelated_xww = intersect(shareid_unrelated_all, shareid_xww) ## 410 unrelated GT individuals who did not fail QC
shareid_related_xww = setdiff(shareid_xww,shareid_unrelated_xww) ## 8187 related GT individuals who did not fail QC

shareid_family2 = read.table(file="../work_Sheng/family_2",header=FALSE)$V2
length(intersect(shareid_related_xww,shareid_family2)) # 8161 related individuals present in both data sets
####shareid_family2: 8173 individuals
####shareid_related_xww: 8187 individuals
####intersection: 8161 individuals
####shareid_family2 and shareid_related_xww are obtained by excluding individuals having genotyping completeness <=96% and estimated inbreeding coefficient>=0.05. 
#### They are slightly different.


#### We use the set "shareid_related_zhw" and "shareid_unrelated_xww" with individuals
####in cohort 2 and with elevated blood glucose level outside the age range 35 and 75
####removed.
## remove the set of individuals in cohort 2 from "shareid_unrelated"
intersect(shareid2_all,shareid_unrelated_xww) # it turns out no one is in cohort 2.
## remove individuals in "shareid_unrelated_xww" whose elevated blood glucose levels are outside our age range 35 and 75 from "shareid_unrelated"
####In this process, I also check whether I could recover Zuohengs phenotypes for 7678 individuals finally retained.

shareid_analysis = c(shareid_related_zhw, shareid_unrelated_xww)
shareid0_analysis = intersect(shareid0_all, shareid_analysis)
shareid1_analysis = intersect(shareid1_all, shareid_analysis)
shareid2_analysis = intersect(shareid2_all, shareid_analysis)
shareid3_analysis = intersect(shareid3_all, shareid_analysis)


#### Determine trait for cohort 0
data_origpheno272 = read.table(file="../work_Sheng/original_pheno_27_2",header=TRUE,
sep="\t",fill=TRUE)
shareid_origpheno272=data_origpheno272$shareid
data_origage = read.table(file="../work_Sheng/original_age",fill=TRUE,header=TRUE,sep="\t")
shareid_origage <- data_origage$shareid
####0: unknown, 1: unaffected, 2: affected, -4: elevated blood glucose level
####outside age range 35 and 75
pheno0_analysis = data.frame(ID = shareid0_analysis, pheno = sapply(shareid0_analysis, function(id){
	if(id %in% shareid_origpheno272){
		bgl <- findInterval(data_origpheno272[shareid_origpheno272 == id, 'age'], c(35,75), rightmost.closed=T)
		if(data_origpheno272[shareid_origpheno272 == id, 'DFDIAB27'] > 0 & bgl == 1) phen <- 2 
		if(data_origpheno272[shareid_origpheno272 == id, 'DFDIAB27'] > 0 & bgl != 1) phen <- -4
		if(data_origpheno272[shareid_origpheno272 == id, 'DFDIAB27'] <= 0 & data_origpheno272[shareid_origpheno272 == id, 'age'] >= 70) phen <- 1
		if(data_origpheno272[shareid_origpheno272 == id, 'DFDIAB27'] <= 0 & data_origpheno272[shareid_origpheno272 == id, 'age'] < 70) phen <- 0
	} else{
		phen <- if(data_origage[shareid_origage == id, 'max_age'] >= 70) 1 else 0
	}
	return(phen)
}))
###check with pheno_zuoheng
sum(sapply( shareid0_related_zhw, function(id){ 
if(data_pheno_zuoheng[shareid_related_zhw == id, 'pheno'] != pheno0_analysis[shareid0_analysis == id, 'pheno']) 1 else 0}))
table(pheno0_analysis$pheno)
shareid0_final = pheno0_analysis[pheno0_analysis$pheno !=-4,'ID']
pheno0 =  pheno0_analysis[pheno0_analysis$pheno !=-4,]


#### Determine trait for cohort 1
data_offpheno72 = read.table(file="../work_Sheng/offspring_pheno_7_2",header=TRUE,
sep="\t",fill=TRUE)
shareid_offpheno72=data_offpheno72$shareid
data_offage = read.table(file="../work_Sheng/offspring_age",fill=TRUE,header=TRUE,sep="\t")
shareid_offage <- data_offage$shareid

pheno1_analysis = data.frame(ID = shareid1_analysis, pheno = sapply(shareid1_analysis, function(id){
	if(id %in% shareid_offpheno72){
		bgl <- findInterval(data_offpheno72[shareid_offpheno72 == id, 'age'], c(35,75), rightmost.closed=T)
		if(data_offpheno72[shareid_offpheno72 == id, 'DFDIAB7'] > 0 & bgl == 1) phen <- 2 
		if(data_offpheno72[shareid_offpheno72 == id, 'DFDIAB7'] > 0 & bgl != 1) phen <- -4
		if(data_offpheno72[shareid_offpheno72 == id, 'DFDIAB7'] <= 0 & data_offpheno72[shareid_offpheno72 == id, 'age'] >= 70) phen <- 1
		if(data_offpheno72[shareid_offpheno72 == id, 'DFDIAB7'] <= 0 & data_offpheno72[shareid_offpheno72 == id, 'age'] < 70) phen <- 0
	} else{
		phen <- if(data_offage[shareid_offage == id, 'max_age'] >= 70) 1 else 0
	}
	return(phen)
}))
###check with pheno_zuoheng
sum(sapply( shareid1_related_zhw, function(id){ 
if(data_pheno_zuoheng[shareid_related_zhw == id, 'pheno'] != pheno1_analysis[shareid1_analysis == id, 'pheno']) 1 else 0}))
table(pheno1_analysis$pheno)
shareid1_final = pheno1_analysis[pheno1_analysis$pheno !=-4,'ID']
pheno1 =  pheno1_analysis[pheno1_analysis$pheno !=-4,]


#### Determine trait for cohort 3
data_gen3pheno = read.table(file="../work_Sheng/gen3_pheno_1_2",header=TRUE,fill=TRUE,sep="\t")
shareid_gen3pheno=data_gen3pheno$shareid
pheno3_analysis = data.frame(ID = shareid3_analysis, pheno = sapply(shareid3_analysis, function(id){
	if(id %in% shareid_gen3pheno){
		bgl <- findInterval(data_gen3pheno[shareid_gen3pheno == id, 'age'], c(35,75), rightmost.closed=T)
		if(data_gen3pheno[shareid_gen3pheno == id, 'DFDIAB1'] > 0 & bgl == 1) phen <- 2 
		if(data_gen3pheno[shareid_gen3pheno == id, 'DFDIAB1'] > 0 & bgl != 1) phen <- -4
		if(data_gen3pheno[shareid_gen3pheno == id, 'DFDIAB1'] <= 0 & data_gen3pheno[shareid_gen3pheno == id, 'age'] >= 70) phen <- 1
		if(data_gen3pheno[shareid_gen3pheno == id, 'DFDIAB1'] <= 0 & data_gen3pheno[shareid_gen3pheno == id, 'age'] < 70) phen <- 0
	}
	return(phen)
}))
###check with pheno_zuoheng
sum(sapply( shareid3_related_zhw, function(id){ 
if(data_pheno_zuoheng[shareid_related_zhw == id, 'pheno'] != pheno3_analysis[shareid3_analysis == id, 'pheno']) 1 else 0}))
table(pheno3_analysis$pheno)
shareid3_final = pheno3_analysis[pheno3_analysis$pheno !=-4,'ID']
pheno3 =  pheno3_analysis[pheno3_analysis$pheno !=-4,]

######################################
#########  *************
####Here is the set of individuals we keep and their phenotype
shareid = c(shareid0_final, shareid1_final,shareid3_final)
pheno = rbind(pheno0,pheno1,pheno3)
shareid0 = intersect(shareid,shareid0_all)
shareid1 = intersect(shareid,shareid1_all)
shareid3 = intersect(shareid,shareid3_all)
length(shareid)
length(shareid0)+length(shareid1)+length(shareid3)

######################################
#########  *************
### Covariate information:
list_c0 = list(
#note that in exam 4, there is no admission blood pressure
#note that height is available in exam 1 and exam 5

## Cohort 0
E7 = rbind(c("MF69","weight_ex1","in pounds"), 
c("MF166","weight_ex2","in pounds"), 
c("MF180","weight_ex3","in pounds"),
c("MF216","weight_ex4","in pounds"),
c("MF292","weight_ex5","in pounds"),
c("MF380","weight_ex6","in pounds"),
c("MF470","weight_ex7","in pounds"),
c("MF67","height_inch_ex1","in full inches"),
c("MF68","height_excess_ex1","in hundredths of inch"),
c("MF290","height_inch_ex5","in full inches"),
c("MF291","height_excess_ex5","in hundredths of inch"),

c("MF56","bp_1_ex1","sys and admission"),
c("MF57","bp_2_ex1","dia and admission"),
c("MF58","bp_3_ex1","sys and examiner 1"),
c("MF59","bp_4_ex1","dia and examiner 1"),
c("MF60","bp_5_ex1","sys and exmainer 2"),
c("MF61","bp_6_ex1","dia and examiner 2"),	
c("MF154","bp_1_ex2","sys and admission"),
c("MF155","bp_2_ex2","dia and admission"),
c("MF156","bp_3_ex2","sys and examiner 1"),
c("MF157","bp_4_ex2","dia and examiner 1"),
c("MF158","bp_5_ex2","sys and exmainer 2"),
c("MF159","bp_6_ex2","dia and examiner 2"),	
c("MF174","bp_1_ex3","sys and admission"),
c("MF175","bp_2_ex3","dia and admission"),
c("MF176","bp_3_ex3","sys and examiner 1"),
c("MF177","bp_4_ex3","dia and examiner 1"),
c("MF178","bp_5_ex3","sys and exmainer 2"),
c("MF179","bp_6_ex3","dia and examiner 2"),	
c("MF264","bp_3_ex4","sys and examiner 1"),
c("MF265","bp_4_ex4","dia and examiner 1"),
c("MF266","bp_5_ex4","sys and exmainer 2"),
c("MF267","bp_6_ex4","dia and examiner 2"),	
c("MF305","bp_1_ex5","sys and admission"),
c("MF306","bp_2_ex5","dia and admission"),
c("MF307","bp_3_ex5","sys and examiner 1"),
c("MF308","bp_4_ex5","dia and examiner 1"),
c("MF309","bp_5_ex5","sys and exmainer 2"),
c("MF310","bp_6_ex5","dia and examiner 2"),		
c("MF384","bp_1_ex6","sys and admission"),
c("MF385","bp_2_ex6","dia and admission"),
c("MF386","bp_3_ex6","sys and examiner 1"),
c("MF387","bp_4_ex6","dia and examiner 1"),
c("MF388","bp_5_ex6","sys and exmainer 2"),
c("MF389","bp_6_ex6","dia and examiner 2"),	
c("MF476","bp_1_ex7","sys and admission"),
c("MF477","bp_2_ex7","dia and admission"),
c("MF478","bp_3_ex7","sys and examiner 1"),
c("MF479","bp_4_ex7","dia and examiner 1"),
c("MF480","bp_5_ex7","sys and exmainer 2"),
c("MF481","bp_6_ex7","dia and examiner 2")		
),

#there is no height measurement
E8 = rbind(c("FA2","weight_ex8","in pounds"),
c("FA7","bp_1_ex8","sys and admission"),
c("FA8","bp_2_ex8","dia and admission"),
c("FA9","bp_3_ex8","sys and examiner 1"),
c("FA10","bp_4_ex8","dia and examiner 1"),
c("FA11","bp_5_ex8","sys and exmainer 2"),
c("FA12","bp_6_ex8","dia and examiner 2")	
),

#there is no height measurement
E9 = rbind(c("FB9","weight_ex9","in pounds"),
c("FB18","bp_1_ex9","sys and admission"),
c("FB19","bp_2_ex9","dia and admission"),
c("FB20","bp_3_ex9","sys and examiner 1"),
c("FB21","bp_4_ex9","dia and examiner 1"),
c("FB22","bp_5_ex9","sys and exmainer 2"),
c("FB23","bp_6_ex9","dia and examiner 2")	
),

E10 = rbind(c("FC6","weight_ex10","in pounds"),
c("FC7","height_ex10","to lower quarter inch"),
c("FC12","bp_1_ex10","sys and admission"),
c("FC13","bp_2_ex10","dia and admission"),
c("FC14","bp_3_ex10","sys and examiner 1"),
c("FC15","bp_4_ex10","dia and examiner 1"),
c("FC16","bp_5_ex10","sys and exmainer 2"),
c("FC17","bp_6_ex10","dia and examiner 2")	
),

E11 = rbind(c("FD12","weight_ex11","in pounds"),
c("FD13","height_ex11","to lower quarter inch"),
c("FD20","bp_1_ex11","sys and admission"),
c("FD21","bp_2_ex11","dia and admission"),
c("FD22","bp_3_ex11","sys and examiner 1"),
c("FD23","bp_4_ex11","dia and examiner 1"),
c("FD24","bp_5_ex11","sys and exmainer 2"),
c("FD25","bp_6_ex11","dia and examiner 2")	
),

E12 = rbind(c("FE67","weight_ex12","in pounds"),
c("FE68","height_ex12","to lower quarter inch"),
c("FE73","bp_1_ex12","sys and admission"),
c("FE74","bp_2_ex12","dia and admission"),
c("FE75","bp_3_ex12","sys and examiner 1"),
c("FE76","bp_4_ex12","dia and examiner 1"),
c("FE77","bp_5_ex12","sys and exmainer 2"),
c("FE78","bp_6_ex12","dia and examiner 2")	
),

E13 = rbind(c("FF63","weight_ex13","in pounds"),
c("FF64","height_ex13","to lower quarter inch"),
c("FF69","bp_1_ex13","sys and admission"),
c("FF70","bp_2_ex13","dia and admission"),
c("FF71","bp_3_ex13","sys and examiner 1"),
c("FF72","bp_4_ex13","dia and examiner 1"),
c("FF73","bp_5_ex13","sys and exmainer 2"),
c("FF74","bp_6_ex13","dia and examiner 2")	
),

E14 = rbind(c("FG62","weight_ex14","in pounds"),
c("FG63","height_ex14","to lower quarter inch"),
c("FG68","bp_1_ex14","sys and admission"),
c("FG69","bp_2_ex14","dia and admission"),
c("FG70","bp_3_ex14","sys and examiner 1"),
c("FG71","bp_4_ex14","dia and examiner 1"),
c("FG72","bp_5_ex14","sys and exmainer 2"),
c("FG73","bp_6_ex14","dia and examiner 2")	
),

E15 = rbind(c("FH62","weight_ex15","in pounds"),
c("FH63","height_ex15","to lower quarter inch"),
c("FH68","bp_1_ex15","sys and admission"),
c("FH69","bp_2_ex15","dia and admission"),
c("FH70","bp_3_ex15","sys and examiner 1"),
c("FH71","bp_4_ex15","dia and examiner 1"),
c("FH72","bp_5_ex15","sys and exmainer 2"),
c("FH73","bp_6_ex15","dia and examiner 2")	
),

#note height has entries
E16 = rbind(c("FI14","weight_ex16","in pounds"),
c("FI15","height_inch_ex16","in inches"),
c("FI16","height_excess_ex16","in hundredth of inch"),
c("FI21","bp_1_ex16","sys and admission"),
c("FI22","bp_2_ex16","dia and admission"),
c("FI23","bp_3_ex16","sys and examiner 1"),
c("FI24","bp_4_ex16","dia and examiner 1"),
c("FI25","bp_5_ex16","sys and exmainer 2"),
c("FI26","bp_6_ex16","dia and examiner 2")	
),

E17 = rbind(c("FJ360","weight_ex17","in pounds"),
c("FJ361","height_ex17","to lower quarter inch"),
c("FJ366","bp_1_ex17","sys and admission"),
c("FJ367","bp_2_ex17","dia and admission"),
c("FJ4","bp_3_ex17","sys and examiner 1"),
c("FJ5","bp_4_ex17","dia and examiner 1"),
c("FJ274","bp_5_ex17","sys and exmainer 2"),
c("FJ275","bp_6_ex17","dia and examiner 2")	
),

E18 = rbind(c("FK47","weight_ex18","in pounds"),
c("FK48","height_ex18","to lower quarter inch"),
c("FK55","bp_1_ex18","sys and admission"),
c("FK56","bp_2_ex18","dia and admission"),
c("FK68","bp_3_ex18","sys and examiner 1"),
c("FK69","bp_4_ex18","dia and examiner 1"),
c("FK369","bp_5_ex18","sys and exmainer 2"),
c("FK370","bp_6_ex18","dia and examiner 2")	
),

#There is no blood pressure measurements from physician2 
E19 = rbind(c("FL009","weight_ex19","in pounds"),
c("FL010","height_ex19","to lower quarter inch"),
c("FL021","bp_1_ex19","sys and admission"),
c("FL022","bp_2_ex19","dia and admission"),
c("FL291","bp_3_ex19","sys and examiner 1"),
c("FL292","bp_4_ex19","dia and examiner 1")
),

#There is no blood pressure measurements from physician2 
E20 = rbind(c("FM19","weight_ex20","in pounds"),
c("FM20","height_ex20","to lower quarter inch"),
c("FM31","bp_1_ex20","sys and admission"),
c("FM32","bp_2_ex20","dia and admission"),
c("FM331","bp_3_ex20","sys and examiner 1"),
c("FM332","bp_4_ex20","dia and examiner 1")	
),


E21 = rbind(c("FN8","weight_ex21","in pounds"),
c("FN9","height_ex21","to lower quarter inch"),
c("FN21","bp_1_ex21","sys and admission"),
c("FN22","bp_2_ex21","dia and admission"),
c("FN310","bp_3_ex21","sys and examiner 1"),
c("FN311","bp_4_ex21","dia and examiner 1"),
c("FN404","bp_5_ex21","sys and exmainer 2"),
c("FN404","bp_6_ex21","dia and examiner 2")	
),

E22 = rbind(c("FO007","weight_ex22","in pounds"),
c("FO008","height_ex22","to lower quarter inch"),
c("FO020","bp_1_ex22","sys and admission"),
c("FO021","bp_2_ex22","dia and admission"),
c("FO329","bp_3_ex22","sys and examiner 1"),
c("FO330","bp_4_ex22","dia and examiner 1"),
c("FO424","bp_5_ex22","sys and exmainer 2"),
c("FO425","bp_6_ex22","dia and examiner 2")	
),

E23 = rbind(c("FP007","weight_ex23","in pounds"),
c("FP008","height_ex23","to lower quarter inch"),
c("FP025","bp_1_ex23","sys and admission"),
c("FP026","bp_2_ex23","dia and admission"),
c("FP296","bp_3_ex23","sys and examiner 1"),
c("FP297","bp_4_ex23","dia and examiner 1"),
c("FP390","bp_5_ex23","sys and exmainer 2"),
c("FP391","bp_6_ex23","dia and examiner 2")	
),

E24 = rbind(c("FQ005","weight_ex24","in pounds"),
c("FQ006","height_ex24","to lower quarter inch"),
c("FQ013","bp_1_ex24","sys examiner 1 first reading"),
c("FQ014","bp_2_ex24","dia examiner 1 first reading"),
c("FQ017","bp_3_ex24","sys examiner 1 second reading"),
c("FQ018","bp_4_ex24","dia examiner 1 second reading"),
c("FQ304","bp_5_ex24","sys and exmainer 2"),
c("FQ305","bp_6_ex24","dia and examiner 2")	
),

E25 = rbind(c("FR005","weight_ex25","in pounds"),
c("FR006","height_ex25","to lower quarter inch"),
c("FR013","bp_1_ex25","sys and admission"),
c("FR014","bp_2_ex25","dia and admission"),
c("FR239","bp_3_ex25","sys and physician first reading"),
c("FR240","bp_4_ex25","dia and physician first reading"),
c("FR380","bp_5_ex25","sys and physician second reading"),
c("FR381","bp_6_ex25","dia and physician second reading")	
),

E26 = rbind(c("FS005","weight_ex26","in pounds"),
c("FS008","height_ex26","to lower quarter inch"),
c("FS015","bp_1_ex26","sys and admission"),
c("FS016","bp_2_ex26","dia and admission"),
c("FS312","bp_3_ex26","sys and physician first reading"),
c("FS313","bp_4_ex26","dia and physician first reading"),
c("FS479","bp_5_ex26","sys and physician second reading"),
c("FS480","bp_6_ex26","dia and physician second reading")
),

E27 = rbind(c("FT004","weight_ex27","in pounds"),
c("FT008","height_ex27","to lower quarter inch"),
c("FT018","bp_1_ex27","sys and admission"),
c("FT019","bp_2_ex27","dia and admission"),
c("FT319","bp_3_ex27","sys and physician first reading"),
c("FT320","bp_4_ex27","dia and physician first reading"),
c("FT483","bp_5_ex27","sys and physician second reading"),
c("FT484","bp_6_ex27","dia and physician second reading")	
)
)

cov0 = cbind(shareid = shareid0)
tot_indiv = nrow(cov0)
for (f in 7:27) {
	file_name = paste("../work_Sheng/clinic_exams/C0_E",f,".txt", sep='')
	data_cov = read.table(file=file_name, sep="\t",fill=TRUE,header=TRUE)
	cov_names = list_c0[[f-6]][,1]
	cov_names_semantic = list_c0[[f-6]][,2]
	len_cov_names = length(cov_names)
	ind_order_in_data <- match(shareid0, data_cov$shareid)
	
	cov_current <- t(sapply(1:tot_indiv, function(i){
		   if (!is.na(ind_order_in_data[i])) {
				unlist(data_cov[ind_order_in_data[i],cov_names])
		   } else rep(NA, len_cov_names)
		   }))
	colnames(cov_current) = cov_names_semantic
	cov0 = cbind(cov0, cov_current)
	print(paste("Exam",f,"is done!"))
}
# Matrix with ID and covariate information
cov0 = as.data.frame(cov0)

ind_all <- match(shareid0, data_shareid_all$shareid)
ind_orig <- match(shareid0,  data_origage$shareid)
pheno0 <- cbind(data_shareid_all$pedno[ind_all], pheno0, data_shareid_all$SEX[ind_all], data_origage$max_age[ind_orig],
apply(data_origage[ind_orig,4:30]-data_origage$max_age[ind_orig], 1, function(x) which(x == 0)))

## Determine height and age information
height_ex1 = cov0$height_inch_ex1+cov0$height_excess_ex1/100
height_ex5 = cov0$height_inch_ex5+cov0$height_excess_ex5/100
height_ex10 = cov0$height_ex10/100
height_ex11 = cov0$height_ex11/100
height_ex12 = cov0$height_ex12/100
height_ex13 = cov0$height_ex13/100
height_ex14 = cov0$height_ex14/100
height_ex15 = cov0$height_ex15/100
height_ex16 = cov0$height_inch_ex16+cov0$height_excess_ex16/100
height_ex17 = cov0$height_ex17/100
height_ex18 = cov0$height_ex18
height_ex19 = cov0$height_ex19
height_ex20 = cov0$height_ex20
height_ex21 = cov0$height_ex21
height_ex22 = cov0$height_ex22
height_ex23 = cov0$height_ex23
height_ex24 = cov0$height_ex24
height_ex25 = cov0$height_ex25
height_ex26 = cov0$height_ex26
height_ex27 = cov0$height_ex27

height = cbind(height_ex1,height_ex5,height_ex10,height_ex11,height_ex12,
height_ex13,height_ex14,height_ex15,height_ex16,height_ex17,
height_ex18,height_ex19,height_ex20,height_ex21,height_ex22,
height_ex23,height_ex24,height_ex25,height_ex26,height_ex27)

thr3=2
# determine rows with extreme values for height
diff3big_set=shareid0[apply(height,1, function(x) sum(diff(sort(x)) > thr3)) > 0]
print(c(thr3,length(diff3big_set)))
# skip modifying following rows (not sure why?)
manual_good = c(291, 485, 634, 665, 775) 

height0_correct = as.data.frame(t(sapply(1:nrow(height), function(i){
	sapply(1:ncol(height), function(j){
		if (!any(abs(height[i,-j]-height[i,j])<=thr3, na.rm=T) & !is.element(i,manual_good)) NA else height[i,j]
	})
})))
colnames(height0_correct) <- colnames(height)
# Update the variables
height_ex2 = height0_correct$height_ex1
height_ex3 = height0_correct$height_ex1
height_ex4 = height0_correct$height_ex5
height_ex6 = height0_correct$height_ex5
height_ex7 = height0_correct$height_ex5
height_ex8 = height0_correct$height_ex10
height_ex9 = height0_correct$height_ex10
# Duplicate some of the columns (ex1, ex5, ex10)
height0_correct = as.data.frame(cbind(height0_correct[,1], height_ex2, 
height_ex3, height_ex4, height0_correct[,2], height_ex6,height_ex7,
height_ex8, height_ex9,height0_correct[,3:20]))
colnames(height0_correct)[c(1,5)] = c("height_ex1","height_ex5")
ncol(height0_correct)

####We work on weight
weight0_correct = as.matrix(cov0[,grepl("weight",colnames(cov0))])

# construct bmi
BMI0_27exams = weight0_correct/height0_correct/height0_correct*703.07
BMI_isna = !is.na(BMI0_27exams)
BMI_lastexam_no <- apply(BMI_isna, 1 , function(x) max(which(x)))

# Get last BMI value on file as well as average BMI for all individuals
BMI_lastexam  <- sapply(1:nrow(BMI0_27exams), function(i) BMI0_27exams[i, BMI_lastexam_no[i]])
BMI_ave = rowMeans(BMI0_27exams, na.rm=T)
# Add BMI variables to data frame
pheno0 = as.data.frame(cbind(pheno0,BMI_lastexam, BMI_lastexam_no,
BMI_ave, 0)) 
colnames(pheno0)=c("pedno","shareid","pheno","sex","lastexam_age",
"lastexam_no","BMI_lastexam","BMI_lastexam_no","BMI_ave","cohortid")

## Cohort 1
list_c1 = list(

E1 = rbind(c("A50","weight_ex1","in pounds"),
c("A51","height_inch_ex1","in full inches"),
c("A52","height_excess_ex1","in hundredths of inch")
),

E2 = rbind(c("B15","weight_ex2","in pounds"),
c("B16","height_inch_ex2","in full inches"),
c("B17","height_excess_ex2","in hundredths of inch")
),

E3 = rbind(c("C416","weight_ex3","in pounds"),
c("C417","height_ex3","to lower quarter inch")
),

E4 = rbind(c("D401","weight_ex4","in pounds"),
c("D402","height_ex4","to lower quarter inch")
),

E5 = rbind(c("E024","weight_ex5","in pounds"),
c("E025","height_ex5","to lower quarter inch")
),

E6 = rbind(c("F007","weight_ex6","in pounds"),
c("F008","height_ex6","to lower quarter inch")
),

E7 = rbind(c("G440","weight_ex7","in pounds"),
c("G441","height_ex7","to lower quarter inch")
)
)

cov1 = cbind(shareid = shareid1)
tot_indiv = nrow(cov1)
for (f in 1:7) {
	file_name = paste("../work_Sheng/clinic_exams/C1_E",f,".txt", sep='')
	data_cov = read.table(file=file_name, sep="\t",fill=TRUE,header=TRUE)
	cov_names = list_c1[[f]][,1]
	cov_names_semantic = list_c1[[f]][,2]
	len_cov_names = length(cov_names)
	ind_order_in_data <- match(shareid1, data_cov$shareid)
	
	cov_current <- t(sapply(1:tot_indiv, function(i){
							if (!is.na(ind_order_in_data[i])) {
							unlist(data_cov[ind_order_in_data[i],cov_names])
							} else rep(NA, len_cov_names)
							}))
	colnames(cov_current) = cov_names_semantic
	cov1 = cbind(cov1, cov_current)
	print(paste("Exam",f,"is done!"))
}
# Matrix with ID and covariate information
cov1 = as.data.frame(cov1)
# Make phenotype+covariate matrix
ind_all <- match(shareid1, data_shareid_all$shareid)
ind_off <- match(shareid1,  data_offage$shareid)
pheno1 <- cbind(data_shareid_all$pedno[ind_all], pheno1, data_shareid_all$SEX[ind_all], data_offage$max_age[ind_off], apply(data_offage[ind_off,4:10]-data_offage$max_age[ind_off], 1, function(x) which(x == 0)))


# Check height information
height_ex1 = cov1$height_inch_ex1+cov1$height_excess_ex1/100
height_ex2 = cov1$height_inch_ex2+cov1$height_excess_ex2/100
height_ex3 = cov1$height_ex3
height_ex4 = cov1$height_ex4
height_ex5 = cov1$height_ex5
height_ex6 = cov1$height_ex6
height_ex7 = cov1$height_ex7
height = cbind(height_ex1,height_ex2, height_ex3, height_ex4, height_ex5,
height_ex6,height_ex7)

thr3=2
# determine rows with extreme values for height
diff3big_set=shareid1[apply(height,1, function(x){ sx <- diff(sort(x)); sum(sx > thr3)>0 | length(sx) == 0})]
print(c(thr3,length(diff3big_set)))
# skip modifying following rows (not sure why?)
manual_good = c(26, 45, 60, 67, 371, 436, 454, 459, 551, 552,
570, 571, 574, 699, 717, 750, 817, 849, 854, 856,
863, 866, 918, 971, 972,1003, 1226, 1388, 1393, 1516,
1625, 1642, 1650, 1829, 1930, 1996, 2002, 2047, 2090, 2092,
2136, 2137, 2144, 2226, 2264, 2392, 2417, 2628, 2682, 2703,
2741, 2827, 2840, 2955, 2970, 3033, 3158, 3174, 3193, 3204, 
3216, 3250, 3252, 3254, 3256, 3259, 3304,3347,3405) 


height1_correct = as.data.frame(t(sapply(1:nrow(height), function(i){
	sapply(1:ncol(height), function(j){
		   if (!any(abs(height[i,-j]-height[i,j])<=thr3, na.rm=T) & !is.element(i,manual_good)) NA else height[i,j]
		   })
})))
colnames(height1_correct) <- colnames(height)

#### weight informaton
weight1_correct = as.matrix(cov1[,grepl("weight",colnames(cov1))])
# construct bmi
BMI1_7exams = weight1_correct/height1_correct/height1_correct*703.07
BMI_isna = !is.na(BMI1_7exams)
BMI_lastexam_no <- apply(BMI_isna, 1 , function(x) max(which(x)))

# Get last BMI value on file as well as average BMI for all individuals
BMI_lastexam  <- sapply(1:nrow(BMI1_7exams), function(i) BMI1_7exams[i, BMI_lastexam_no[i]])
BMI_ave = rowMeans(BMI1_7exams, na.rm=T)
# Add BMI variables to data frame
pheno1 = as.data.frame(cbind(pheno1,BMI_lastexam, BMI_lastexam_no,
BMI_ave, 1)) 
colnames(pheno1)=c("pedno","shareid","pheno","sex","lastexam_age",
"lastexam_no","BMI_lastexam","BMI_lastexam_no","BMI_ave","cohortid")

## Cohort 3
data_cov = read.csv("../work_Sheng/clinic_exams/C3_E1.csv",fill=TRUE,header=TRUE)
ind_order_in_data <- match(shareid3, data_cov$shareid)
weight3 = data_cov[ind_order_in_data, "G3A444"]
height3 = data_cov[ind_order_in_data, "G3A446"]
bmi3 = data_cov[ind_order_in_data, "G3A707"]

####four individuals have missing height and weight
BMI3 = weight3/height3/height3*703.07
which(is.na(bmi3))
which(is.na(BMI3))
# Compare BMI from data to that computed "by hand"
max(abs(BMI3-bmi3),na.rm=TRUE) #not quite the same -- rounding??

# Make phenotype+covariate matrix
ind_all <- match(shareid3, data_shareid_all$shareid)
ind_off <- match(shareid3,  data_gen3pheno$shareid)
pheno3 <- cbind(data_shareid_all$pedno[ind_all], pheno3, data_shareid_all$SEX[ind_all], data_gen3pheno$age_exam[ind_off], 1)
# Add BMI variables to data frame
pheno3 = as.data.frame(cbind(pheno3,BMI3, 1, BMI3,3))
colnames(pheno3)=c("pedno","shareid","pheno","sex","lastexam_age",
"lastexam_no","BMI_lastexam","BMI_lastexam_no","BMI_ave","cohortid")

# Note that in the old version, only BMI average across exams was recorded in pheno (eg pheno_sheng_v1)
pheno = rbind(pheno0,pheno1,pheno3)
pheno$pedno[is.na(pheno$pedno)] <- 0
dim(pheno)
### Will use sex and average BMI as covariates in null model
write.table(pheno,"phenocov.txt",sep="\t", row.names=FALSE,
col.names=TRUE)






