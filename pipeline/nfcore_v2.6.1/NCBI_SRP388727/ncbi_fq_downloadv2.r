library(dplyr)
#FQ_ID<-c("ERR1018186","ERR1018189","ERR1018190","ERR1018192","ERR1018193","ERR1018194","ERR1018199","ERR1018202","ERR1018205","ERR1018206","ERR1018210","ERR1018213","ERR1018215","ERR1018216","ERR1018217","ERR1018218","ERR1018219","ERR1018222","ERR1018223","ERR1018224","ERR1018225","ERR1018226","ERR1018227","ERR1018228","ERR1018229","ERR1018239","ERR1018245","ERR1018246","ERR1018248","ERR1018249","ERR1018262","ERR1018273","ERR1018274","ERR1018276","ERR1018287")
FQ_ID=read.csv("SRR_Acc_List.txt",col.names = FALSE)
FQ_ID=FQ_ID$FALSE. %>% unique()
print("inputID.length:")
print(length(FQ_ID))

for(i in FQ_ID){
  all_file<-list.files()
  #print(i)
  if(grepl(i,all_file) %>% any()!=TRUE){
    print(i)
    print("download FQ")
    cmd1<-sprintf(' /home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s  | grep "ftp.sra" | sed "s/\\"url\\"://g" | sed -n 1p |xargs  curl -O',i)
    cmd2<-sprintf(' /home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s  | grep "ftp.sra" | sed "s/\\"url\\"://g" | sed -n 2p | xargs  curl -O',i)
    print(cmd1)
    print(cmd2)
    system(cmd1)
    system(cmd2)


    print("download md5")
    system("touch all.md5")
    cmd3<-sprintf('/home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s | jq -r ".[] | \\"\\(.md5) \\(.filename)\\"" >> all.md5',i)
    # cmd3<-sprintf(' /home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s  | jq -r ".[] | "\\(.filename) \\(.md5)" " >> md5.id.txt ',i)
    print(cmd3)
    system(cmd3)
    
    # system('''ffq --ftp $ID | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs echo| cut -d " " -f1 | xargs curl -O''')
    # system('''ffq --ftp $ID | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs echo| cut -d " " -f2 | xargs curl -O''')
  }
  
}
#
#获取md5方法1
#for(i in FQ_ID){
#  # i<-"SRR14610569"
#  all_file<-list.files()
#  if(grepl(i,all_file) %>% any() !=FALSE ){
#    print(i)
#    print("download md5")
#    
#    # get md5
#    system("touch all.md5")
#    cmd3<-sprintf('/home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s | jq -r ".[] | \\"\\(.md5) \\(.filename)\\"" >> all.md5',i)
#    # cmd3<-sprintf(' /home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s  | jq -r ".[] | "\\(.filename) \\(.md5)" " >> md5.id.txt ',i)
#    print(cmd3)
#    system(cmd3)
#    
#    # # get filename
#    # system("touch md5.filename.txt")
#    # cmd5<-sprintf(' /home/zhiyu/data/software/anaconda3/bin/ffq --ftp %s  | jq -r ".[] | .filename" >> md5.filename.txt ',i)
#    # print(cmd5)
#    # system(cmd5)    
#    
#  }
#  
#}


#paste md5.id.txt md5.filename.txt > md5.all.txt

#获取md5方法2
##pip3 install jtbl
system("cat all.md5  | parallel --pipe -N1  md5sum -c >md5_result")

