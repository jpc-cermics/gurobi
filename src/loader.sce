// this loader is not generated 
// we have added the global=%t option 
// because some other dynamic libraries 
// may want to use symbols defined internally here.

libgurobi_path=file('join',['.','libgurobi'+%shext]);
addinter(libgurobi_path,'libgurobi',global=%t);
clear libgurobi_path;
