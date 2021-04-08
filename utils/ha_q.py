import os
from pbs_jobs import submit


haq_path = os.path.join(os.path.expanduser('~'), '.HaQ')
with open(haq_path, 'r') as handle:
    cmds = handle.readlines()
job_ran = True
while job_ran and len(cmds)>0:
    cmd = cmds.pop(0)
    job_ran = submit(cmd, queue=False)
if len(cmds) > 0:
    cmds = [cmd] + cmds
with open(haq_path, 'w') as handle:
        handle.writelines(cmds)

