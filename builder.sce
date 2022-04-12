// build scicos toolbox

libname='nspgurobi'

// generate Path.incl file with relative path to nsp
// ------------------------------------------------

ilib_path_incl(relative=%f)

// compile shared library in src.
// ------------------------------

if c_link('libgurobi_Interf') then
  printf('please do not rebuild a shared library while it is linked\n')
  printf('use ulink to unlink first\n');
else 
  printf('building shared library\n')
  chdir('src');
  // we need to chdir to execute a builder.
  ok=exec('builder.sce',errcatch=%t);
  if ~ok then 
    x_message('Compilation of source file failed\n");
  end
  chdir("../");
end 

// macros 
//--------

// add_lib('macros',compile=%t);

// [4] man 
// chdir('man')
// exec('builder.sce') 
// chdir(dir);




