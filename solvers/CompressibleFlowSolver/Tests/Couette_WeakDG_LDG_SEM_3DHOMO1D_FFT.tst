<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000397706</value>
            <value variable="rhou" tolerance="1e-12">48.1295</value>
            <value variable="rhov" tolerance="1e-12">0.145836</value>
             <value variable="rhow" tolerance="1e-12">8.79677e-06</value>
            <value variable="E" tolerance="1e-12">17519.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00139611</value>
            <value variable="rhou" tolerance="1e-12">83.3516</value>
            <value variable="rhov" tolerance="1e-12">0.505196</value>
             <value variable="rhow" tolerance="1e-12">3.17716e-05</value>
            <value variable="E" tolerance="1e-12">18953</value>
        </metric>
    </metrics>
</test>


