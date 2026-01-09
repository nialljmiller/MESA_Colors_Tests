import os
import shutil
import subprocess

mesa_dir = os.environ["MESA_DIR"]

# cp synthetic.f90 $MESA_DIR/colors/private/
shutil.copy(
    "synthetic.f90",
    os.path.join(mesa_dir, "colors/private/synthetic.f90"),
)

# cd $MESA_DIR
# ./clean
subprocess.run(["./clean"], cwd=mesa_dir, check=True)

# ./install
subprocess.run(["./install"], cwd=mesa_dir, check=True)

