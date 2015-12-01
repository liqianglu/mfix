import thread
import time
from twisted.internet import reactor, protocol

import mfix

def run_mfix():
    mfix.mfix()

class Echo(protocol.Protocol):
    """This is just about the simplest possible protocol"""

    def dataReceived(self, data):
        "As soon as any data is received, write it back."
        self.transport.write(data)


def main():
    """This runs the protocol on port 8000"""
    factory = protocol.ServerFactory()
    factory.protocol = Echo
    reactor.listenTCP(8000,factory)
    reactor.run()
    thread.start_new_thread(run_mfix, ())

# this only runs if the module was *not* imported
if __name__ == '__main__':
    main()
