from logging import *

# Removes all log handlers

#try:
handlers = Logger.manager.loggerDict#.keys()

for name, hdlr in handlers.iteritems():
        Logger.removeHandler(getLogger(), hdlr)
        
print 'All log handlers removed'
#except:
#    print 'Could not remove handlers.'

log = getLogger()

x = list(log.handlers)
for i in x:
    log.removeHandler(i)
    i.flush()
    i.close()
    
print Logger.manager.loggerDict.keys()

