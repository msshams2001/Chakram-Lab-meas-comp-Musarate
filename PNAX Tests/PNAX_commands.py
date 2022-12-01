class PNAX():
    def __init__():


    def set_pulse_output(self, ch, n, value, name):
        """
        Turns the pulse output ON and OFF.
        Parameters
        
        <ch> Any existing channel number; if unspecified, value is set to 1.

        <n> Pulse generator number. Choose from 0 to 4. 0 is the generator that pulses the ADC.

        <value> Boolean

                ON (or 1) - turns pulse output ON.

                OFF (or 0) - turns pulse output OFF.

        <name>   Optional. String name of the pulse generator.

                Required for use with external pulse generators.

        Use SENSe:PULSe:CAT? to return the names of configured pulse generators.

    If specified, <n> is ignored.

    If unspecified, <n> is required for internal pulse generators.

    Examples

    SENS:PULS1 1

    SENS:PULS 1, "My81110"

    Query Syntax

    SENSe<ch>:PULSe[:STATe]? [<name>]

    Return Type

    Boolean

    Default

    OFF
        """
        
        query = "SENS%d:PULS%d[:STAT] %d,%s" %(ch,n,value,name)
        return self.write(query)


    def set_pulse_delay(self, ch, n, value, name = None):
        '''
        <ch>Any existing channel number; if unspecified, value is set to 1.


        <n> Internal pulse generator number. Choose from 0 to 4. 0 is the generator that pulses the ADC.

        <value> Delay value in seconds. Choose a value from about 33ns to about 70 seconds.

        <name> Optional. String name of the pulse generator.

        Required for use with external pulse generators.

        Use SENSe:PULSe:CAT? to return the names of configured pulse generators.

        If specified, <n> is ignored.

        If unspecified, <n> is required for internal pulse generators.


        SENSe<ch>:PULSe<n>:DELay <value>[,<name>]
        '''

        if ch is None:
            ch = 1

        query = "SENS%d:PULS%d:DEL %f, '%s'" % (ch,n,value.name)
        return self.write(query)
