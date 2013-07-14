<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC, supersonic Navier-Stokes equations parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Perturbation_NS_M15_square_CBC_par.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> Perturbation_NS_M15_square_CBC_par.xml</file>
        <file description="Restart File"> Perturbation_NS_M15_square_CBC_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">6.78852e-06</value>
            <value variable="rhou" tolerance="1e-12">0.0025743</value>
            <value variable="rhov" tolerance="1e-12">0.00155051</value>
            <value variable="E" tolerance="1e-12">17373.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.000231713</value>
            <value variable="rhou" tolerance="1e-12">0.211816</value>
            <value variable="rhov" tolerance="1e-12">0.0502772</value>
            <value variable="E" tolerance="1e-12">69537.9</value>
        </metric>
    </metrics>
</test>


