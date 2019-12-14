from tkinter import Label, Entry, Listbox,Scrollbar, END, VERTICAL, N, S, StringVar
 
class NoteFrm:
    frm_base=None
    
    label_list=[]
    entry_list=[]
    data_list=[]
    
    totn=2
    row=False

    def __init__(self,frm_base,init_notes=[],init_data=[0,0],totn=2,row=False,entryWidth=8,labelWidth=8):
        self.frm_base=frm_base
        self.totn=len(init_data)
        self.row=row
        self.entryWidth=entryWidth
        self.labelWidth=labelWidth
        self.loadlabel(frm_base,init_notes)
        self.loadentries(frm_base,init_data)

        
        
    def loadlabel(self,frm_base,init_notes):
        self.label_list=[]
        for i in range(self.totn):
            label=Label(frm_base,text=init_notes[i])
            self.label_list.append(label)
        
        if self.row:
            for i in range(self.totn):
                self.label_list[i].grid(row=0,column=2*i)
                #self.label_list[i+1].grid(row=0,column=3*i+2)
        else:
            for i in range(self.totn):
                self.label_list[i].grid(row=i,column=0)
                #self.label_list[2*i+1].grid(row=i,column=2)
    
        
    def loadentries(self,frm_base,init_data):
        self.entry_list=[]
        for i in range(self.totn):
            entry=Entry(frm_base,width=self.entryWidth)
            entry.insert(0,init_data[i])
            self.entry_list.append(entry)
        
        if self.row:
            for i in range(self.totn):
                self.entry_list[i].grid(row=0,column=2*i+1)

        else:
            for i in range(self.totn):
                self.entry_list[i].grid(row=i,column=1)
    
    def setValue(self,data=[1,600]):
        #print data
        for i in range(self.totn):
            self.entry_list[i].delete(0,END)
            self.entry_list[i].insert(0,data[i])
        
    def report(self,strrep=False):
        self.data_list=[]
        if strrep:
            for ele in self.entry_list:
                value=ele.get()
                #value=value if value else '0'
                self.data_list.append(value)
        else:
            for ele in self.entry_list:
                value=ele.get()
                self.data_list.append(float(value))
        return self.data_list

class InfoList:

    def __init__(self,baseFrame, selectEvent, dataFormat, initialInfo=[],width=30,height=5):
        self.infoData=initialInfo
        self.dataFormat=dataFormat
        self.infoList=Listbox(baseFrame,width=width,height=height,listvariable=StringVar(value=[ele[0] for ele in self.infoData]))
        self.infoList.bind('<<ListboxSelect>>',selectEvent)
        self.infoList.grid(row=0,column=0)

        self.scrollbar=Scrollbar(baseFrame,orient=VERTICAL, command=self.infoList.yview)
        self.scrollbar.grid(row=0,column=1,sticky=(N,S))
        self.infoList['yscrollcommand']=self.scrollbar.set

        self.updateInfo(self.infoData)

    def updateInfo(self,infoData=[]):
        self.infoList.delete(0,len(self.infoData)-1)
        self.infoData=infoData
        for iinfo, info in enumerate(self.infoData):
            info[0]=iinfo
            self.infoList.insert(iinfo,self.dataFormat(info))
        for i in range(0,len(self.infoData),2):
            self.infoList.itemconfigure(i, background='#f0f0ff')

    def report(self):
        idxs = self.infoList.curselection()
        if len(idxs)==1:
            return self.infoData[idxs[0]]
        return []