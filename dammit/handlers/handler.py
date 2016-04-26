from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

class TaskHandler(TaskLoader):

    def __init__(self, args, doit_db, files=None, **doit_config_kwds):
        super(TaskHandler, self).__init__()

        if files is None:
            self.files = {}
        if type(files) is not dict:
            raise TypeError('files must be of type dict')

        else:
            self.files = files

        self.tasks = {}
        self.doit_config = dict(dep_file=doit_db, **doit_config_kwds)
        self.doit_dep_mgr = D


    def register_task(name, task, files=None):
        if files is None:
            self.files = {}
        if type(files) is not dict:
            raise TypeError('files must be of type dict')
        
        self.tasks[name] = task
        self.files.update()

    
    def load_tasks(self, cmd, opt_values, pos_args):
        return self.tasks.values()

    def run(self, doit_args=None):
        if doit_args is None:
            doit_args = ['run']
        return DoitMain(self()).run(doit_args)
