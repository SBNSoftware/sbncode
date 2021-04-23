
class test_msg:

    level=0 # level [0,1,2,3] = [debug,info,warning,error]

    def __init__(self):
        pass

    @classmethod
    def debug(cls,msg=''):
        if int(cls.level)<=0:
            print('\033[94m[DEBUG]  \033[00m',msg)

    @classmethod
    def info(cls,msg=''):
        if int(cls.level)<=1:
            print('\033[92m[INFO]   \033[00m',msg)
    
    @classmethod
    def warning(cls,msg=''):
        if int(cls.level)<=2:
            print('\033[95m[WARNING]\033[00m',msg)
    
    @classmethod
    def error(cls,msg=''):
        if int(cls.level)<=3:
            print('\033[91m[ERROR]  \033[00m',msg)

debug   = test_msg.debug
info    = test_msg.info
warning = test_msg.warning
error   = test_msg.error

