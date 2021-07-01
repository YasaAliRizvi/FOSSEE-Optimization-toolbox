/*
This file prints and saves the outputs of all the problems in the CUTEst folder.

Written by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

funcprot(0);
//List of all CUTEst problems (including rosenbrock which is not a CUTEst problem)
problem_list2 = list("broydn3d", "cosine", "curly10", "dtoc1l", "hadamard", "hager1", "hager4", "liswet10", "mancino", "msqrtb", "penalty1", "power", "reading2", "rosenbrock", "sipow3", "tfi2", "ubh1");
k=1;
temp =1;
i=1;
while i<18
   try
    k=i;
    temp = k;
    printf("\n\nProblem %d  =  ", i);
    p = problem_list2(temp);
    disp(p);
    exec(get_absolute_file_path("run_all.sce")+"\"+ p+".sce", -1);
    printf("fval  =  ");
    disp(fval);
    printf("exitflag  =  ");
    disp(exitflag);
    printf("output  =  ");
    disp(output);
    //To save each problem outputs
    //save(p+".sod","fval","output","exitflag","x");
    
    //To save all problem outputs in one file: (load the file to get the outputs of all the problems)
    //Ex: format for fval of hager1 problem is hager1_fval (similar format for other outputs and problems)
    execstr(p+"_fval = fval");
    execstr(p+"_output = output");
    execstr(p+"_exitflag = exitflag");
    execstr(p+"_x = x");
    save("outputs_all.sod",p+"_fval",p+"_output",p+"_exitflag",p+"_x");
    
    i = temp+1;
    
    catch
        [errmsg, tmp, nline, func] = lasterror()
        msg = "%s: error on line #%d: ""%s""\n"
        msg = msprintf(msg, func, nline, errmsg)
        error(msg);
        i=temp+1;
    end 
    
end
