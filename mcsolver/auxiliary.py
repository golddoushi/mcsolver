'''
Created on 2019 7 6

@author: Andrew
'''
import numpy as np
import os
import shutil

def getMSD(input,overavg=True,bestChoice=True,verbose=False):
    '''
    '''
    #check the shape
    shape=input.shape
    if(len(shape)==1):
        # the case of one line of data
        dim=int(shape[0])
        avg=1.*np.sum(input)/dim
        sd=0.
        diffList=[]
        for ele in input:
            diff=(ele-avg)*(ele-avg)
            diffList.append(diff)
            sd+=diff
        rmsd=np.sqrt(sd)/dim/abs(avg) if overavg else np.sqrt(sd)/dim
        if verbose:
            print('dimension is: %d'%dim)
            print('average is: %.6f'%avg)
            print('relative mean-squared-diff is: %.6f'%rmsd)
        if bestChoice:
            sortList=quicksort(diffList)
            if verbose:
                print('the best chois is:',input[sortList[0]],' with derviation: ', diffList[sortList[0]])
            return rmsd,sortList[0]
        
    elif(len(shape)==2):
        # the case of 2x2 matrix
        ncol=int(shape[1])
        rmsdList=[]
        bestChoiceList=[]
        for icol in range(ncol):
            rmsd,bestChoiceID=getMSD(input[:,icol],bestChoice=True,verbose=verbose)
            rmsdList.append(rmsd)
            bestChoiceList.append(bestChoiceID)
        if bestChoice:
            return rmsdList,bestChoiceList
        else:
            return rmsdList
    else:
        print('cannot handle data dimension >2')
        exit()

def getCombination(inputList=[],ncomb=1):
    '''
    This function can generate combinations from extracted elements in
    inputList, e.g., if we want to select two elements from [1,2,3], 
    they maybe either [1,2], [1,3] or [2,3], thus we can use this facility as
    comb=getCombination(inputList=[1,2,3],ncomb=2)
    then we can get a list of: [ [1,2], [1,3], [2,3] ]
    
    params definition:
    inputList: the pool of elements
    ncomb    : num. of elements in combination
    '''
    #print 'conduct function'
    #print 'input:', inputList
    ntot=len(inputList)
    if(ncomb==1):
        result=[]
        for ele in inputList:
            result.append([ele])
        return result
    else:
        if(ncomb>ntot):
            print('Error: num. of combinations larger than input!')
            exit()
        else:
            result=[]
            for ihead in range(0,ntot-ncomb+1):
                headele=inputList[ihead]
                #print 'head element:',headele
                resList=list(inputList[ihead+1:ntot])
                #print 'reside list:',resList
                resCombList=getCombination(resList,ncomb=ncomb-1)
                for ele in resCombList:
                    ele.insert(0,headele)
                    #print 'ele',ele
                    result.append(list(ele))
                    #print 'ele',ele
                #result.append(result_head)
            return result

def doesTheTwoListHaveSameNumber(list1,list2):
    for n1 in list1:
        for n2 in list2:
            if n1==n2:
                return True
    return False

def findTheFirstSameNumberAmongTwoList(list1,list2): # if no same element presented, return -1
    for n1 in list1:
        for n2 in list2:
            if n1==n2:
                return n1
    return -1

def improveTheMatrixRankToThree(mat:np.ndarray):
    if np.linalg.matrix_rank(mat)==3: return mat

    if np.linalg.matrix_rank(mat)==2:
        # try to fill one of the dimension
        for idim in range(3):
            mat_trial=np.array(list(mat))
            mat_trial[idim,idim]=1
            if np.linalg.matrix_rank(mat_trial)==3: return mat_trial
        raise('Fail to improve the matrix rank')
    
    if np.linalg.matrix_rank(mat)==1:
        for idim in range(2):
            for jdim in range(idim+1,3):
                mat_trial=np.array(list(mat))
                mat_trial[idim,idim]=1
                mat_trial[jdim,jdim]=1
                if np.linalg.matrix_rank(mat_trial)==3: return mat_trial
        raise('Fail to improve the matrix rank')

class Node:
    def __init__(self,value=0,maxValue=9,leftNode=False):
        self.value=value
        self.maxValue=maxValue
        self.leftNode=leftNode
    def addone(self):
        self.value+=1
        if self.value>self.maxValue:
            # carry site
            if self.leftNode and self.leftNode.addone():
                self.value=self.leftNode.value+1
            else:
                return False
        if self.leftNode:
            self.leftNode.maxValue=self.value-1
        return True
    
def getCombinationLoop(input_list=[],ncomb=1,maxComb=-1):
    totInput=len(input_list)
    index_list=[i for i in range(totInput)]
    # create linked list
    rootNode=Node(value=0,maxValue=totInput-1)
    for i in range(1,ncomb):
        node_i=Node(value=i,maxValue=totInput-1,leftNode=rootNode)
        rootNode=node_i
    # generate combination numbers in one loop
    def getTree(remoteNode):
        tree=[remoteNode.value]
        while(remoteNode.leftNode):
            tree.insert(0,remoteNode.leftNode.value)
            remoteNode=remoteNode.leftNode
        return tree
    comb_list=[[input_list[i] for i in getTree(rootNode)]]
    ncount=0
    while(rootNode.addone()):
        comb_list.append([input_list[i] for i in getTree(rootNode)])
        #print(getTree(rootNode))
        ncount+=1
        if maxComb>0 and ncount>=maxComb:
            break
    return comb_list


def quicksort(input_list=[]):
    '''
    Sort the input list and return a sort index list
    e.g. we want to sort input_list = [1, 3, 2, 5, 0]
    then the function will return the sorted INDEX:
        output = [4, 0, 2, 1, 3]
        its first element is 4, since the lowest num. 0
        is at 4th position in original input_list, and so on
    
    WARNING: the input_list will also be changed forever. 
    '''
    def __firstSort(input_list,sort_list,low,high):
        l=low;h=high
        split=input_list[low]
        
        while(l<h):
            while(input_list[h]>split):
                h-=1
                
            if(l<h):
                input_list[l]=input_list[h]
                input_list[h]=split
                
                sort_tmp=sort_list[l]
                sort_list[l]=sort_list[h]
                sort_list[h]=sort_tmp
                
                l+=1
            
            while(input_list[l]<split):
                l+=1
            
            if l<h :
                input_list[h]=input_list[l]
                input_list[l]=split
                
                sort_tmp=sort_list[l]
                sort_list[l]=sort_list[h]
                sort_list[h]=sort_tmp
                
                h-=1
        
        if l>low :
            __firstSort(input_list,sort_list,low,h-1)
            
        if h<high:
            __firstSort(input_list,sort_list,l+1,high)
    if len(input_list)==1:
        return [0]    
    if len(input_list)==0:
        print('ERROR: try to sort empty list')
        exit()
        return []
    sort_list=[]
    high=len(input_list)-1
    for i in range(0,high+1):
        sort_list.append(i)
    __firstSort(input_list,sort_list,0,high)            
    return sort_list
        
def createDir(path):
    try:
        os.mkdir(path)
    except OSError as err:
        print(err,'--> dir %s exist, skip'%path)
        #print('WARNING: %s exist'%path)

def ls(path):
    try:
        return os.listdir(path)
    except OSError as err:
        print(err, '--> cannot find %s, stop'%path)
        exit()

def cp(file_src, file_des):
    try:
        shutil.copy2(file_src,file_des)
    except IOError as err:
        print(err, '--> copy file %s to %s failed, stop'%(file_src,file_des))
        exit()
    #except FileNotFoundError as err:
    #    print(err, '--> cannot find %s, stop'%file_src)
    #    exit()
    #except shutil.SameFileError as err:
    #    print(err, '--> %s and %s are the same file, skip'%(file_src,file_des))

def rm(path,quiet=False):
    try:
        os.remove(path)
    except OSError as err:
        if not quiet: print(err, '--> remove file %s failed, skip'%path)
