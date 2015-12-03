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

        # make sure that data is exactly 512 characters in size
        data = data[:512]
        data = list(data)+['']*(512-len(data))

        # wait for semaphore to be free
        while pymfix.main.in_semaphore[0]:
            pass
        pymfix.main.in_semaphore[0] = True

        pymfix.main.in_buffer = data

        # wait for mfix to reply
        pymfix.main.out_semaphore[0] = True
        while pymfix.main.out_semaphore[0]:
            pass

        # print response from mfix
        print(''.join(pymfix.main.out_buffer))
        self.transport.write(''.join(pymfix.main.out_buffer))
        pymfix.main.in_semaphore[0] = False

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
