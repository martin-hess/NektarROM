<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady FRDG advection GAUSS_LAGRANGE, P=10, homogeneous Dirichlet bcs, deformed elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_ISO_deformed_GAUSS_LAGRANGE_3x3.xml</parameters>
    <files>
        <file description="Session File">Advection2D_ISO_deformed_GAUSS_LAGRANGE_3x3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.0248746 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.0951207 </value>
        </metric>
    </metrics>
</test>
