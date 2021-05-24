#!/usr/bin/env python3
"""
Tools to generate LSF commands
"""

def bsub(command, stdout=None, stderr=None, ram=None, group=None,
         name=None, dep=None, cwd=None, gpu=0, gpu_exclusive=False,
         queue=None, project=None, hosts=None):
    """Generate LSF submission command string"""

    # Construct resource usage string
    rusage = {}
    if ram:
        rusage['mem'] = ram

    rusage_str = ', '.join([f'{key}={value}' for key, value in rusage.items()])

    # Construct gpu string
    if gpu == 1 and not gpu_exclusive:
        gpu = '-'
    elif gpu > 0:
        gpu = f'"num={gpu}:j_exclusive={"yes" if gpu_exclusive else "no"}"'

    # Construct command
    command = ['bsub',
               f'-J "{name}"' if name else '',
               f'-g "{group}"' if group else '',
               f'-P "{project}"' if project else '',
               f'-q "{queue}"' if queue else '',
               f'-m "{hosts}"' if hosts else '',
               f'-w "{dep}"' if dep else '',
               f'-gpu {gpu}' if gpu else '',
               f'-M {ram}' if ram else '',
               f'-R "rusage[{rusage_str}]"' if rusage else '',
               f'-o {stdout}' if stdout else '',
               f'-e {stderr}' if stderr else '',
               f'-cwd "{cwd}"' if cwd else '',
               f'"{command}"']

    command = [x for x in command if x]

    return ' '.join(command)
