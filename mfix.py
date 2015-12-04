import thread
import time
from twisted.internet import reactor, protocol

import pymfix

def run_mfix():
    pymfix.main.mfix()

class Echo(protocol.Protocol):
    """This is just about the simplest possible protocol"""

    def dataReceived(self, data):
        "Called with data sent from client."
        pymfix.main.request_pending = 1

        # wait for mfix to stop
        while 1 != pymfix.main.mfix_waiting:
            pass

        # change any variable...
        # pymfix.run.time  = 1000
        pymfix.run.tstop = pymfix.run.tstop/2
        # pymfix.run.dt    = 10 # has no effect?

        response = 'the pressure is %s \n tstop is %s' % (str(pymfix.fldvar.p_g[:100]),pymfix.run.tstop)
        # print(response)
        self.transport.write(response)

        # resume mfix
        pymfix.main.request_pending = 0

def main():
    """This runs the protocol on port 8000"""
    thread.start_new_thread(run_mfix, ())
    factory = protocol.ServerFactory()
    factory.protocol = Echo
    reactor.listenTCP(8000,factory)
    reactor.run()

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
