import os
from setuptools import setup, find_packages

f = open(os.path.join("darpanx","get_dir.py"), "w")
f.write("import os"+ '\n')
f.write("dir=r'"+os.path.join(os.getcwd(),"darpanx'")+ '\n')
f.write("def get_dir():"+ '\n')
f.write(" return dir"+ '\n')
f.write("def get_wrap():"+ '\n')
f.write(" return os.path.join(dir,'darpanx_fit.py')"+ '\n')
f.close()

setup(name='DarpanX',
      version='0.3',
      description='DarpanX: A Python Package for Modeling X-ray Reflectivity of Multilayer Mirrors',
      author='Biswajit Mondal',
      author_email='biswajitm@prl.res.in',
      url='https://github.com/biswajitmb/DarpanX/tree/objOriented/objectOriented/python',
      packages=find_packages()
     )

print("%% DarpanX_message : Don't delete the directory '"+os.getcwd()+"'")

